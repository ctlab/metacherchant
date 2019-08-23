package io;

import org.apache.log4j.Logger;
import ru.ifmo.genetics.dna.Dna;
import ru.ifmo.genetics.io.ReadersUtils;
import ru.ifmo.genetics.io.sources.NamedSource;
import ru.ifmo.genetics.structures.map.BigLong2ShortHashMap;
import ru.ifmo.genetics.utils.NumUtils;
import ru.ifmo.genetics.utils.tool.ExecutionFailedException;
import ru.ifmo.genetics.utils.tool.Tool;
import utils.HashFunction;
import utils.PolynomialHash;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ConcurrentMap;
import java.util.concurrent.CountDownLatch;

import static utils.StringUtils.normalizeDna;

/**
 * Created by -- on 23.08.2019.
 */
public class LargeKmerLoader {
    static final int READS_WORK_RANGE_SIZE = 1 << 15;   // 32 K

    public static HashFunction hash;

    static class ReadsLoadWorker extends ReadsWorker {
        ReadsLoadWorker(Map<String, Integer> subgraph, int k, int minDnaLen, HashFunction hasher, BigLong2ShortHashMap graph) {
            this.subgraph = subgraph;
            this.graph = graph;
            this.k = k;
            this.minDnaLen = minDnaLen;
            this.hasher = hasher;
        }

        final BigLong2ShortHashMap graph;
        final Map<String, Integer> subgraph;
        final int k;
        final int minDnaLen;
        final HashFunction hasher;

        @Override
        public void process(List<Dna> reads) {
            for (Dna dna : reads) {
                if (dna.length() >= minDnaLen) {
                    for (int i = 0; i + k <= dna.length(); i++) {
                        String key = dna.substring(i, i + k).toString();
                        subgraph.put(normalizeDna(key), (int) graph.get(hasher.hash(normalizeDna(key))));
                    }
                }
            }
        }
    }

    public static Map<String, Integer> loadReads(File[] files, int k, int minSeqLen, int availableProcessors,
                                                 Logger logger, BigLong2ShortHashMap graph, HashFunction hasher)
            throws ExecutionFailedException {
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 12, true);
        ConcurrentMap<String, Integer> subgraph = new ConcurrentHashMap<>();
        hash = hasher != null ? hasher : new PolynomialHash();

        LargeKmerLoader.ReadsLoadWorker[] workers = new LargeKmerLoader.ReadsLoadWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new LargeKmerLoader.ReadsLoadWorker(subgraph, k, minSeqLen, hash, graph);
        }

        run(files, workers, hm, logger);

        return subgraph;
    }

    public static void run(File[] files, ReadsWorker[] workers, BigLong2ShortHashMap hmForMonitoring, Logger logger)
            throws ExecutionFailedException {
        for (File file : files) {
            Tool.info(logger, "Loading file " + file.getName() + "...");

            NamedSource<Dna> reader = null;
            try {
                reader = ReadersUtils.readDnaLazyTrunc(file, null);
            } catch (IOException e) {
                throw new ExecutionFailedException("Failed to read from file " + file.getPath());
            }

            ReadsDispatcher dispatcher = new ReadsDispatcher(reader, READS_WORK_RANGE_SIZE, hmForMonitoring);
            CountDownLatch latch = new CountDownLatch(workers.length);

            for (int i = 0; i < workers.length; ++i) {
                workers[i].setDispatcher(dispatcher);
                workers[i].setLatch(latch);
                new Thread(workers[i]).start();
            }

            try {
                latch.await();
            } catch (InterruptedException e) {
                Tool.warn(logger, "Main thread interrupted");
                for (ReadsWorker worker : workers) {
                    worker.interrupt();
                }
                throw new ExecutionFailedException("Thread was interrupted", e);
            }
            Tool.info(logger, NumUtils.groupDigits(dispatcher.reads) + " reads added");
        }
    }
}

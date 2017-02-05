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

import java.io.*;
import java.util.List;
import java.util.concurrent.CountDownLatch;

public class  LargeKIOUtils {

    static final int READS_WORK_RANGE_SIZE = 1 << 15;   // 32 K

    public static HashFunction hash;

    static class ReadsLoadWorker extends ReadsWorker {
        ReadsLoadWorker(BigLong2ShortHashMap hm, int k, int minDnaLen, HashFunction hasher) {
            this.hm = hm;
            this.k = k;
            this.minDnaLen = minDnaLen;
            this.hasher = hasher;
        }

        final BigLong2ShortHashMap hm;
        final int k;
        final int minDnaLen;
        final HashFunction hasher;
        int totalSeq = 0, goodSeq = 0;
        long totalLen = 0, goodLen = 0;

        @Override
        public void process(List<Dna> reads) {
            for (Dna dna : reads) {
                totalSeq++;
                totalLen += dna.length();

                if (dna.length() >= minDnaLen) {
                    for (int i = 0; i + k <= dna.length(); i++) {
                        long hash = hasher.hash(dna, i, i + k);
                        hm.addAndBound(hash, (short) 1);
                    }
                    goodSeq++;
                    goodLen += dna.length();
                }
            }
        }
    }

    public static BigLong2ShortHashMap loadReads(File[] files, int k, int minSeqLen,
                                                 int availableProcessors, Logger logger)
            throws ExecutionFailedException {
        BigLong2ShortHashMap hm = new BigLong2ShortHashMap(
                (int) (Math.log(availableProcessors) / Math.log(2)) + 4, 12, true);

        ReadsLoadWorker[] workers = new ReadsLoadWorker[availableProcessors];
        for (int i = 0; i < workers.length; ++i) {
            workers[i] = new ReadsLoadWorker(hm, k, minSeqLen, hash != null ? hash : new PolynomialHash());
        }

        run(files, workers, hm, logger);

        // calculating statistics...
        int totalSeq = 0, goodSeq = 0;
        long totalLen = 0, goodLen = 0;
        for (ReadsLoadWorker worker : workers) {
            totalSeq += worker.totalSeq;
            goodSeq += worker.goodSeq;
            totalLen += worker.totalLen;
            goodLen += worker.goodLen;
        }
        Tool.debug(logger,
                "Good/Total sequences count = " + NumUtils.groupDigits(goodSeq) + "/" + NumUtils.groupDigits(totalSeq)
                + " (" + String.format("%.1f", goodSeq * 100.0 / totalSeq) + "%)");
        Tool.debug(logger,
                "Good/Total sequences length = " + NumUtils.groupDigits(goodLen) + "/" + NumUtils.groupDigits(totalLen)
                        + " (" + String.format("%.1f", goodLen * 100.0 / totalLen) + "%)");
        logger.debug("k-mers HM size = " + NumUtils.groupDigits(hm.size()));

        return hm;
    }

    public static void run(File[] files, ReadsWorker[] workers, BigLong2ShortHashMap hmForMonitoring, Logger logger)
            throws ExecutionFailedException {
        for (File file : files) {
            Tool.info(logger, "Loading file " + file.getName() + "...");

            NamedSource<Dna> reader = ReadersUtils.readDnaLazy(file);

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

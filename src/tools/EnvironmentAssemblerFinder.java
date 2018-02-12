package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Map;

import ru.ifmo.genetics.utils.tool.Parameter;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.FileParameterBuilder;
import ru.ifmo.genetics.utils.tool.inputParameterBuilder.StringParameterBuilder;

public class EnvironmentAssemblerFinder {
    public final Parameter<String> assembler = new Parameter(new StringParameterBuilder("assembler")
            .withDescription("assembler wich you want to use")
            .withDefaultValue("None")
            .create());

    public final Parameter<String> assemblerPath = new Parameter(new StringParameterBuilder("assemblerpath")
            .withDescription("path of the assembler")
            .withDefaultValue("")
            .create());

    public final Parameter<File> outputDir = new Parameter(new FileParameterBuilder("output")
            .mandatory()
            .withShortOpt("o")
            .withDescription("output directory")
            .create());


    private static String getOutputPrefix(Map<String, String> dictionaryArgs) {
        String outputPrefix = "";
        if (dictionaryArgs.containsKey("--output")) outputPrefix = dictionaryArgs.get("--output");
        if (dictionaryArgs.containsKey("-o")) outputPrefix = dictionaryArgs.get("-o");
        outputPrefix += "/";
        return outputPrefix;
    }

    public static Map<String, String> runAssembler(Map<String, String> dictionaryArgs) {
        String assemb = dictionaryArgs.get("--assembler"); //assembler.get();
        String outputPrefix = getOutputPrefix(dictionaryArgs);
        Map<String, String> dicAfterAssemblerArgs = new HashMap<String, String>();// = dictionaryArgs;
        dicAfterAssemblerArgs.put("--seq", dictionaryArgs.get("--seq"));
        dicAfterAssemblerArgs.put("--coverage", "0");
        dicAfterAssemblerArgs.put("--maxkmers", dictionaryArgs.get("--maxkmers"));
        dicAfterAssemblerArgs.put("-o", outputPrefix + "/result");


        if (assemb.compareTo("spades") == 0) {
            try {
                String[] command = {"python", dictionaryArgs.get("--assemblerpath") + "/spades.py", "--12",
                        outputPrefix + "cutReads.fastq", "-o", outputPrefix + "out_spades"};
                ProcessBuilder procBuilder = new ProcessBuilder(command);
                procBuilder.redirectErrorStream(true);
                Process process = procBuilder.start();
                InputStream stdout = process.getInputStream();
                InputStreamReader isrStdout = new InputStreamReader(stdout);
                BufferedReader brStdout = new BufferedReader(isrStdout);
                String line = null;
                while ((line = brStdout.readLine()) != null) {
                    System.out.println(line);
                }

                dicAfterAssemblerArgs.put("--reads", outputPrefix + "out_spades/" + "contigs.fasta");
                // or better take assembly_graph.gfa?
                dicAfterAssemblerArgs.put("--k", "55");
            } catch (IOException e) {
                // TODO Auto-generated catch block
                e.printStackTrace();
            }
            return dicAfterAssemblerArgs;
        }
        if (assemb.compareTo("megahit") == 0) {
            try {
                String[] commandA = {dictionaryArgs.get("--assemblerpath") + "/megahit", "--12",
                        outputPrefix + "cutReads.fastq", "-o", outputPrefix + "out_megahit"};
                String[] commandR = {"mv", outputPrefix + "out_megahit/" + "final.contigs.fa",
                        outputPrefix + "out_megahit/" + "final.contigs.fasta"};
                ProcessBuilder procBuilder = new ProcessBuilder(commandA);
                ProcessBuilder procBuilderR = new ProcessBuilder(commandR);

                procBuilder.redirectErrorStream(true);
                Process process = procBuilder.start();
                InputStream stdout = process.getInputStream();
                InputStreamReader isrStdout = new InputStreamReader(stdout);
                BufferedReader brStdout = new BufferedReader(isrStdout);

                String line = null;
                while ((line = brStdout.readLine()) != null) {
                    System.out.println(line);
                }

                procBuilderR.redirectErrorStream(true);
                Process processB = procBuilderR.start();
                InputStream stdoutB = processB.getInputStream();
                InputStreamReader isrStdoutB = new InputStreamReader(stdoutB);
                BufferedReader brStdoutB = new BufferedReader(isrStdoutB);
                String lineB = null;
                while ((lineB = brStdoutB.readLine()) != null) {
                    System.out.println(lineB);
                }

                dicAfterAssemblerArgs.put("--reads", outputPrefix + "out_megahit/" + "final.contigs.fasta");
                dicAfterAssemblerArgs.put("--k", "55");
            } catch (IOException e) {
                e.printStackTrace();
            }
            return dicAfterAssemblerArgs;
        }
        dicAfterAssemblerArgs.put("--k", "55");
        return dicAfterAssemblerArgs;
    }

    public static void main(String[] args) {
        new EnvironmentFinderMain().mainImpl(args);


        Map<String, String> dictionaryArgs = new HashMap<String, String>();
        for (int i = 0; i < args.length - 1; i += 2) {
            System.out.println(args[i]);
            dictionaryArgs.put(args[i], args[i + 1]);
        }

        dictionaryArgs = runAssembler(dictionaryArgs);

        String[] argsAfterAssembler = new String[dictionaryArgs.keySet().size() * 2];

        for (int i = 0; i < dictionaryArgs.keySet().size(); i++) {
            String key = dictionaryArgs.keySet().toArray()[i].toString();
            argsAfterAssembler[i * 2] = key;
            argsAfterAssembler[i * 2 + 1] = dictionaryArgs.get(key);
        }

        new EnvironmentFinderMain().mainImpl(argsAfterAssembler);

    }
}
//--k 3 --reads reads.fasta --seq genes.fasta -o ress --assembler megahit --assemblerpath /home/mariya/megahit/megahit --maxkmers 1000 --coverage 0 --procfiltration 100

package io;

import ru.ifmo.genetics.dna.DnaQ;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

public class RichFastaReader {

    private BufferedReader br;

    public RichFastaReader(File f) throws FileNotFoundException {
        br = new BufferedReader(new FileReader(f));
    }

    public RichFastaReader(String s) throws FileNotFoundException {
        this(new File(s));
    }

    private List<DnaQ> dnas;
    private List<String> comments;

    public List<DnaQ> getDnas() {
        if (dnas == null) {
            readAll();
        }
        return dnas;
    }

    public List<String> getComments() {
        if (dnas == null) {
            readAll();
        }
        return comments;
    }

    private void readAll() {
        dnas = new ArrayList<DnaQ>();
        comments = new ArrayList<String>();

        try {
            boolean lastComment = true;
            StringBuilder curComment = new StringBuilder();
            StringBuilder curDna = new StringBuilder();
            while (br.ready()) {
                String line = br.readLine();
                if (line.startsWith(">") || line.startsWith(";")) {
                    if (!lastComment) {
                        dnas.add(new DnaQ(curDna.toString(), 0));
                        curDna.setLength(0);
                        curComment.setLength(0);
                    }
                    line = line.substring(1);
                    curComment.append(line);
                    lastComment = true;
                } else {
                    if (lastComment) {
                        comments.add(curComment.toString());
                        curDna.setLength(0);
                        curComment.setLength(0);
                    }
                    curDna.append(line);
                    lastComment = false;
                }
            }
            if (curComment.length() > 0) {
                comments.add(curComment.toString());
            }
            if (curDna.length() > 0) {
                dnas.add(new DnaQ(curDna.toString(), 0));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}

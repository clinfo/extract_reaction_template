package src.main;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.standardizer.Standardizer;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import com.chemaxon.mapper.AutoMapper;
import com.chemaxon.mapper.Mapper;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;

public class extractReactionTemplate {
    @Option(name = "-i", aliases = {"--input"}, metaVar = "input", required = true, usage = "Path to input reaction file.")
    private static String INPUT_FILE;

    @Option(name = "-o", aliases = {"--output"}, metaVar = "output", required = true, usage = "Path to output.")
    private static String OUTPUT_FILE;

    @Option(name = "-f", aliases = {"--format"}, metaVar = "format", required = true,
            usage = "Output file format. (ex: smiles, smarts, sdf, etc.)")
    private static String FILE_FORMAT;

    public static void main(String[] args) {
        extractReactionTemplate ert = new extractReactionTemplate();
        CmdLineParser parser = new CmdLineParser(ert);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.out.println("USAGE: (All arguments is required.)");
            System.out.println();
            parser.printUsage(System.out);
            return;
        }
        try(BufferedReader br = Files.newBufferedReader(Paths.get(INPUT_FILE), StandardCharsets.UTF_8)){
            String line;
            String[] splitWords;
            List<InputStream> molStreams = new ArrayList<>();
            List<Integer> patentNumberList = new ArrayList<>();
            ByteArrayInputStream is;
            int i = 0;
            while (((line = br.readLine()) != null)) {
                splitWords = line.split("\t");
                if (splitWords.length < 1) {
                    continue;
                }
                String smiles = splitWords[0];
                i++;
                //String patentNo = splitWords[1];
                is = new ByteArrayInputStream(smiles.getBytes(StandardCharsets.UTF_8));
                molStreams.add(is);
                patentNumberList.add(i);
            }
            writeListToOutput(molStreams, patentNumberList, OUTPUT_FILE);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static void writeListToOutput(List<InputStream> molStreams, List<Integer> idList, String save_path) throws IOException{
        try(BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(save_path), StandardCharsets.UTF_8))) {
            Standardizer std = new Standardizer(new File("./strip_salts.xml"));
            MolImporter mi;
            Molecule mol;
            writer.write("id,product,reaction,reaction_center,reaction_1neighbor");
            writer.newLine();
            for (int i = 0; i < molStreams.size(); i++) {
                try {
                    mi = new MolImporter(molStreams.get(i));
                    mol = mi.read();
                    mi.close();
                    RxnMolecule reaction = RxnMolecule.getReaction(mol);
                    Utils.stripSalts(reaction, std);
                    if (reaction.getProductCount() != 1) {
                        continue;
                    }
                    AutoMapper.unmap(reaction);
                    Molecule product = reaction.getProduct(0);
                    AutoMapper mapper = new AutoMapper();
                    mapper.setMappingStyle(Mapper.MappingStyle.CHANGING);
                    mapper.map(reaction);
                    RxnMolecule reactionCloneForCore = reaction.clone();
                    RxnMolecule reactionCloneFor1Neighbor = reaction.clone();
                    Core.getReactionTemplate(reactionCloneForCore, "0");
                    Core.getReactionTemplate(reactionCloneFor1Neighbor, "1");
                    if (reactionCloneForCore.getProduct(0).isEmpty() | reactionCloneFor1Neighbor.getProduct(0).isEmpty()) {
                        continue;
                    }
                    AutoMapper.unmap(product);
                    AutoMapper.unmap(reaction);
                    AutoMapper.unmap(reactionCloneForCore);
                    AutoMapper.unmap(reactionCloneFor1Neighbor);
                    reactionCloneForCore.removeEmptyComponents();
                    reactionCloneFor1Neighbor.removeEmptyComponents();
                    writer.write(idList.get(i) + "," +
                            MolExporter.exportToFormat(product, FILE_FORMAT) + "," +
                            MolExporter.exportToFormat(reaction, FILE_FORMAT) + "," +
                            MolExporter.exportToFormat(reactionCloneForCore, FILE_FORMAT) + "," +
                            MolExporter.exportToFormat(reactionCloneFor1Neighbor, FILE_FORMAT));
                    writer.newLine();
                } catch (ArrayIndexOutOfBoundsException | IllegalArgumentException | MolFormatException e) {
                    e.printStackTrace();
                }
            }
            System.out.println("[SAVE] " + save_path);
        }
    }

}


package src.main;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import com.chemaxon.mapper.AutoMapper;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.BufferedWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

public class AtomMappingOnReactionTemplate {
    @Option(name = "-i", aliases = {"--input"}, metaVar = "input", required = true, usage = "Path to input reaction template file")
    private static String INPUT_FILE;

    public static void main(String[] args) {
        AtomMappingOnReactionTemplate amrt = new AtomMappingOnReactionTemplate();
        CmdLineParser parser = new CmdLineParser(amrt);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            parser.printUsage(System.out);
            System.exit(1);
        }
        List<String> output = new ArrayList<>();
        try(MolImporter importer = new MolImporter(INPUT_FILE)) {
            Molecule mappedReaction;
            // AutoMapper configuration
            AutoMapper.Options options = new AutoMapper.Options();
            options.setMappingStyle(AutoMapper.MappingStyle.COMPLETE);
            //
            Iterator<Molecule> mappedMoleculeIterator = AutoMapper.iterator(importer.iterator(), options);
            while (mappedMoleculeIterator.hasNext()) {
                mappedReaction = RxnMolecule.getReaction(mappedMoleculeIterator.next());
                output.add(MolExporter.exportToFormat(mappedReaction, "smarts"));
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        int lastDot = INPUT_FILE.lastIndexOf(".");
        String output_filename = INPUT_FILE.substring(0, lastDot) + "_mapped" + INPUT_FILE.substring(lastDot);
        Path output_path = Paths.get(output_filename);
        try (BufferedWriter writer = Files.newBufferedWriter(output_path)) {
            for (String line : output) {
                writer.append(line);
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

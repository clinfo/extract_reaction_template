package src.main;

import chemaxon.formats.MolImporter;
import chemaxon.standardizer.Standardizer;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import com.chemaxon.mapper.AutoMapper;
import com.chemaxon.mapper.Mapper;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;

import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.List;


class Utils {
    static Molecule getMolFromInputStream(InputStream is) {
        Molecule mol = null;
        try (MolImporter mi = new MolImporter(is)) {
            mol = mi.read();
        } catch (IOException e) {
            e.printStackTrace();
            System.exit(1);
        }
        return mol;
    }

    static void mapReaction(RxnMolecule rxn, Mapper.MappingStyle ms) {
        AutoMapper mapper = new AutoMapper();
        mapper.setMappingStyle(ms);
        mapper.map(rxn);
    }

    static void parseArgument(Object instance, String[] args) {
        CmdLineParser parser = new CmdLineParser(instance);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            System.out.println("USAGE: (All arguments is required.)");
            System.out.println();
            parser.printUsage(System.out);
            System.exit(1);
        }
    }

    static Boolean checkMwOfProductsInReaction(RxnMolecule reaction) throws IOException {
        double threshold = 1000;
        List<Boolean> booleanList = new ArrayList<>();
        for (int i = 0; i < reaction.getProductCount(); i++) {
            booleanList.add(reaction.getProduct(i).getMass() > threshold);
        }
        return booleanList.contains(true);
    }

    static void stripSalts(RxnMolecule reaction, Standardizer std) {
        int j = reaction.getComponentCount(RxnMolecule.AGENTS);
        for (int i = 0; i < j; i++) {
            reaction.removeComponent(RxnMolecule.AGENTS, 0);
        }
        for (int i = 0; i < reaction.getReactantCount(); i++) {
            Molecule reactant = reaction.getReactant(i);
            std.standardize(reactant);
        }
        for (int i = 0; i < reaction.getProductCount(); i++) {
            Molecule product = reaction.getProduct(i);
            std.standardize(product);
        }
        reaction.removeEmptyComponents();
    }
}

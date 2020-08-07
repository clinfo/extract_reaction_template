package src.main;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolFormatException;
import chemaxon.formats.MolImporter;
import chemaxon.marvin.io.formats.mdl.MolExport;
import chemaxon.reaction.Reaction;
import chemaxon.reaction.ReactionException;
import chemaxon.reaction.Reactor;
import chemaxon.sss.SearchConstants;
import chemaxon.sss.search.SearchOptions;
import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import chemaxon.util.iterator.IteratorFactory;
import com.chemaxon.mapper.AutoMapper;
import com.chemaxon.mapper.Mapper;
import org.apache.commons.collections4.map.PassiveExpiringMap;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;

public class filterValidReactionTemplate {
    @Option(name = "-i", aliases = {"--input"}, metaVar = "input", required = true,
            usage = "Path to input file")
    private static String INPUT_FILE;

    public static void main(String[] args) {
        filterValidReactionTemplate fvrt = new filterValidReactionTemplate();
        CmdLineParser parser = new CmdLineParser(fvrt);
        try {
            parser.parseArgument(args);
        } catch (CmdLineException e) {
            parser.printUsage(System.out);
            System.exit(1);
        }

        String[] header = null;
        List<String> outputList = new LinkedList<>();
        try (BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(INPUT_FILE)))) {
            String line;
            RxnMolecule rxnTemplate;
            String sRxnComp;
            int count = 0;
            String[] data;
            while ((line = br.readLine()) != null) {
                if (count == 0) {
                    header = line.split("\t");
                } else {
                    data = line.split("\t");  // data: (product, reaction template, year)
                    rxnTemplate = complementReactionTemplate(data[1]);
                    if (! isValidReaction(data[0], rxnTemplate)) {
                        continue;
                    }
                    sRxnComp = MolExporter.exportToFormat(rxnTemplate, "smarts:u");
                    outputList.add(data[0] + "\t" + sRxnComp + "\t" + data[2]);
                }
                count += 1;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        //
        int lastDotIdx = INPUT_FILE.lastIndexOf(".");
        String outputFilename = INPUT_FILE.substring(0, lastDotIdx) + "_filtered" + INPUT_FILE.substring(lastDotIdx);
        Path outputPath = Paths.get(outputFilename);
        try (BufferedWriter writer = Files.newBufferedWriter(outputPath)) {
            writer.append(String.join("\t", header));
            writer.newLine();
            for (String line : outputList) {
                writer.append(line);
                writer.newLine();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    private static Boolean isValidReaction(String product, RxnMolecule rxnTemplate) {
        SearchOptions searchOptions = new SearchOptions(SearchConstants.DEFAULT_SEARCHTYPE);
        searchOptions.setStereoModel(SearchOptions.STEREO_MODEL_LOCAL);
        searchOptions.setStereoSearchType(SearchOptions.STEREO_EXACT);
        Molecule[] reactantsForward;
        Molecule[] reactantsReverse;
        List<Molecule[]> reactantsList = new ArrayList<>();
        String proUniq;
        String genProduct;
        Reactor reactor;
        try {
            Molecule[] pro = {MolImporter.importMol(product)};
            try {
                proUniq = MolExporter.exportToFormat(pro[0], "smiles:u");
            } catch (IOException e) {
                return false;
            }
            // reverse
            reactor = new Reactor();
            reactor.setReverse(true);
            reactor.setReaction(rxnTemplate);
            reactor.setReactants(pro);
            reactor.setSearchOptions(searchOptions.toString());
            while ((reactantsReverse = reactor.react()) != null) {
                reactantsList.add(reactantsReverse);
            }
            if (reactantsList.contains(null)) {
                return false;
            }
            // forward
            for (Molecule[] reactants : reactantsList) {
                reactor = new Reactor();
                reactor.setReverse(false);
                reactor.setReaction(rxnTemplate);
                reactor.setReactants(reactants);
                reactor.setSearchOptions(searchOptions.toString());
                while ((reactantsForward = reactor.react()) != null) {
                    if (reactantsForward.length != 1) {
                        continue;
                    }
                    try {
                        genProduct = MolExporter.exportToFormat(reactantsForward[0], "smiles:u");
                    } catch (IOException e) {
                        continue;
                    }
                    if (genProduct.equals(proUniq)) {
                        return true;
                    }
                }
            }
           return false;
        } catch (ReactionException | MolFormatException | ArrayIndexOutOfBoundsException e) {
            return false;
        }
    }

    private static RxnMolecule complementReactionTemplate(String sRxn) throws MolFormatException {
        IteratorFactory.AtomIterator iterator;
        Set<Integer> aMapsPro = new HashSet<>();
        Set<Integer> aMapsRcts = new HashSet<>();

        RxnMolecule rxn = RxnMolecule.getReaction(MolImporter.importMol(sRxn));
        atomMappingOnReaction(rxn);
        if (rxn.getProductCount() != 1) {
            System.out.println("[WARNING] A reaction template has only one product template.");
            return null;
        }
        // for product
        Molecule pro = rxn.getProduct(0);
        iterator = pro.getAtomIterator();
        while (iterator.hasNext()) {
            MolAtom atom = iterator.next();
            aMapsPro.add(atom.getAtomMap());
        }
        // for reactants
        for (Molecule rct : rxn.getReactants()) {
            iterator = rct.getAtomIterator();
            while (iterator.hasNext()) {
                MolAtom atom = iterator.next();
                aMapsRcts.add(atom.getAtomMap());
            }
        }
        //
        Set<Integer> aMapsNotInRcts = new HashSet<>(aMapsPro);
        aMapsNotInRcts.retainAll(aMapsRcts);
        //
        Molecule addMol = pro.clone();
        removeAtomsWithAtomMaps(addMol, aMapsNotInRcts);
        //AutoMapper.unmap(rxn);
        try {
            String sMol = MolExporter.exportToFormat(addMol, "smarts");
            if (sMol.contains(".")) {
                return rxn;
            }
            int atomCount = addMol.getAtomCount();
            Set<String> prohibitSymbols = new HashSet<>();
            Collections.addAll(prohibitSymbols, "*", null);
            if (atomCount == 0 || (atomCount == 1 &&prohibitSymbols.contains(addMol.getAtom(0).getSymbol()))) {
                return rxn;
            }
        } catch (IOException e) {
            return rxn;
        }
        rxn.addComponent(addMol, RxnMolecule.REACTANTS);
        //AutoMapper.unmap(rxn);
        return rxn;
    }

    private static void removeAtomsWithAtomMaps(Molecule mol, Set<Integer> aMaps) {
        IteratorFactory.AtomIterator iterator = mol.getAtomIterator();
        MolAtom atom;
        while (iterator.hasNext()) {
            atom = iterator.next();
            if (aMaps.contains(atom.getAtomMap())) {
                iterator.remove();
            }
        }
    }

    private static void atomMappingOnReaction(RxnMolecule rxn) {
        AutoMapper mapper = new AutoMapper();
        mapper.setMappingStyle(Mapper.MappingStyle.COMPLETE);
        mapper.map(rxn);
    }
}

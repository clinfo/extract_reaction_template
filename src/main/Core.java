package src.main;

import chemaxon.struc.MolAtom;
import chemaxon.struc.Molecule;
import chemaxon.struc.RxnMolecule;
import chemaxon.util.iterator.IteratorFactory;
import org.apache.commons.lang3.BooleanUtils;

public class Core {
    static void getReactionTemplate(RxnMolecule reaction, String neighborNum) throws IllegalArgumentException {
        // Process reactant(s) included in a reaction
        for (int i = 0; i < reaction.getReactantCount(); i++) {
            Molecule reactant = reaction.getReactant(i);
            IteratorFactory.AtomIterator iterator = setPropertyOnMolecule(reactant, neighborNum);
            removeAtomWithNoProperty(iterator);
        }
        // Process product(s) included in a reaction
        for (int i = 0; i < reaction.getProductCount(); ++i) {
            Molecule product = reaction.getProduct(i);
            IteratorFactory.AtomIterator iterator = setPropertyOnMolecule(product, neighborNum);
            removeAtomWithNoProperty(iterator);
        }
    }

    private static void removeAtomWithNoProperty(IteratorFactory.AtomIterator iterator) {
        while (iterator.hasNext()) {
            MolAtom atom = iterator.next();
            Boolean reserveAtom = (Boolean) atom.getProperty("reserve");
            Boolean reserveNAtom = (Boolean) atom.getProperty("nReserve");
            if (BooleanUtils.isNotTrue(reserveAtom) & BooleanUtils.isNotTrue(reserveNAtom)) {
                iterator.remove();
            }
        }
    }

    private static IteratorFactory.AtomIterator setPropertyOnMolecule(Molecule mol, String neighborNum) {
        IteratorFactory.AtomIterator iterator;
        switch (neighborNum) {
            case "0": {
                iterator = mol.getAtomIterator();
                setPropertyOnMappedAtom(iterator);
                return mol.getAtomIterator();
            }
            case "1": {
                iterator = mol.getAtomIterator();
                setPropertyOnMappedAtom(iterator);
                iterator = mol.getAtomIterator();
                setPropertyOnNeighborAtom(iterator, mol);
                return mol.getAtomIterator();
            }
            default: {
                System.out.println("[ERROR] check the neighbor atom number is valid");
                System.exit(1);
                return null;  // Avoid compile errors
            }
        }
    }

    private static void setPropertyOnMappedAtom(IteratorFactory.AtomIterator iterator) {
        while (iterator.hasNext()) {
            MolAtom atom = iterator.next();
            int mapNum = atom.getAtomMap();
            if (mapNum != 0) {
                atom.putProperty("reserve", true);
            }
        }
    }

    private static void setPropertyOnNeighborAtom(IteratorFactory.AtomIterator iterator, Molecule reaction) {
        Boolean reserveAtom;
        Boolean reserveNAtom;
        while (iterator.hasNext()) {
            MolAtom atom = iterator.next();
            reserveAtom = (Boolean) atom.getProperty("reserve");
            // For only processing reaction center atoms which atom-mapped by AutoMapper.
            if (BooleanUtils.isTrue(reserveAtom)) {
                // set "reserve" to the reaction center atoms and "nReserve" to the first neighbor atoms.
                IteratorFactory.AtomNeighbourIterator nIterator = new IteratorFactory(reaction).createAtomNeighbourIterator(atom);
                while (nIterator.hasNext()) {
                    MolAtom nAtom = nIterator.next();
                    reserveAtom = (Boolean) nAtom.getProperty("reserve");
                    reserveNAtom = (Boolean) nAtom.getProperty("nReserve");
                    // for preventing overwrite on reaction center atoms / first neighbor atoms
                    if (BooleanUtils.isTrue(reserveAtom) & BooleanUtils.isTrue(reserveNAtom)) {
                        continue;
                    }
                    nAtom.putProperty("nReserve", true);
                }
            }
        }
    }
}

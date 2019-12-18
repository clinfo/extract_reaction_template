package src.main;

import chemaxon.formats.MolExporter;
import chemaxon.formats.MolImporter;
import chemaxon.standardizer.Standardizer;
import chemaxon.struc.Molecule;

import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

import chemaxon.struc.RxnMolecule;
import com.chemaxon.mapper.AutoMapper;
import com.chemaxon.mapper.Mapper.MappingStyle;
import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Option;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import javax.xml.xpath.XPath;
import javax.xml.xpath.XPathConstants;
import javax.xml.xpath.XPathExpressionException;
import javax.xml.xpath.XPathFactory;
import java.io.*;
import java.util.*;
import java.util.regex.Pattern;

import static src.main.Core.getReactionTemplate;
import static src.main.Utils.*;


public class extractReactionTemplateFromReaxysData {
    @Option(name = "-t", aliases = {"--process_type"}, metaVar = "process type", required = true,
            usage = "multiThread or singleThread")
    private static String PROCESS_TYPE = null;

    @Option(name = "-f", aliases = {"--format"}, metaVar = "format", required = true,
            usage = "Output mol format. (ex: smiles, smarts, sdf, etc.)")
    private static String FORMAT = null;

    @Option(name = "-d", aliases = {"--debug"}, metaVar = "debug",
            usage = "Debug mode")
    private static Boolean USE_DEBUG_MODE = false;

    @Option(name = "-m", aliases = {"--process_multi_dirs"}, metaVar = "process multi dirs",
            usage = "Process multiple directories")
    private static Boolean PROC_MULTI_DIRS = false;

    @Option(name = "-p", aliases = {"--process_dir_path"}, metaVar = "process directory", required = true,
            usage = "Absolute path to the directory for processing")
    private static String PROC_DIR_PATH = null;


    public static void main(String[] args) {
        String cacheSize = String.format("%d", 1024 * 1024);
        System.setProperty("chemaxon.automapper.AutoMapper.cacheSize", cacheSize);
        extractReactionTemplateFromReaxysData p = new extractReactionTemplateFromReaxysData();
        parseArgument(p, args);
        File[] dirs = new File[1];
        dirs[0] = new File(PROC_DIR_PATH);
        List<Job> jobs = new ArrayList<>();
        File[] xml_files = getFileNamesFromDir(dirs[0].toString());
        for (File xml_file : xml_files) {
            jobs.add(new Job(xml_file, args));
        }
        if (PROCESS_TYPE.equals("multiThread")) {
            System.out.println("[INFO] multi thread process");
            jobs.parallelStream().forEach(Job::doJob);
        }
        if (PROCESS_TYPE.equals("singleThread")) {
            System.out.println("[INFO] single thread process");
            jobs.forEach(Job::doJob);
        }
    }

    static class Job {
        private File file;
        private String[] args;
        Job(File file, String[] args) {
            Job.this.file = file;
            Job.this.args = args;
        }
        void doJob() {
            try {
                Document document = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(file);
                convertXmlFilesToCsv("//RY.STR", document, "rn", file, args);
                Path from = Paths.get(file.getPath());
                Path to = Paths.get(file.getParent() + "/xml/" + file.getName());
                Files.move(from, to);
            } catch (ParserConfigurationException | SAXException | IOException e) {
                e.printStackTrace();
            }
        }
    }

    private static void convertXmlFilesToCsv(String tagLocation, Document doc, String attrName, File f, String[] args) {
        try {
            XPath xpath = XPathFactory.newInstance().newXPath();
            NodeList nodes = (NodeList) xpath.evaluate(tagLocation, doc, XPathConstants.NODESET);

            List<InputStream> molStreams = new ArrayList<>();
            List<Integer> idList = new ArrayList<>();
            InputStream is;
            for (int i = 0; i < nodes.getLength(); i++) {
                Element element = (Element) nodes.item(i);
                String text = element.getTextContent();
                String attr = element.getAttribute(attrName);
                if (StringUtils.isBlank(attr) || StringUtils.isBlank(text)) {
                    continue;
                }
                is = new ByteArrayInputStream(text.getBytes(StandardCharsets.UTF_8));
                molStreams.add(is);
                idList.add(Integer.parseInt(attr));
            }
            writeListToOutput(molStreams, idList, f);
        } catch (XPathExpressionException e) {
            e.printStackTrace();
        }
    }

    private static File[] getFileNamesFromDir(String homePath) {
        File dir = new File(homePath);
        FilenameFilter filter = (file, name) -> name.endsWith(".xml");
        File[] childFiles = dir.listFiles(filter);
        if (USE_DEBUG_MODE) {
            if (childFiles == null) {
                System.err.println(String.format("[ERROR] %s has no files whose name end with 'xml'", homePath));
                System.exit(1);
            } else {
                for (File f : childFiles) {
                    System.out.println("[PROC] " + f);
                }
            }
        }
        return childFiles;
    }

    private static void writeListToOutput(List<InputStream> molStreams, List<Integer> idList, File xml_file) {
        String[] absPathRmExt;
        if (PROC_MULTI_DIRS) {
            absPathRmExt = xml_file.getPath().split(Pattern.quote("."));
        } else {
            absPathRmExt = xml_file.getAbsolutePath().split(Pattern.quote("."));
        }
        String[] parentDirPathList = xml_file.getParent().split("/");
        String max_pub_year = parentDirPathList[parentDirPathList.length -1];
        String save_path = absPathRmExt[0] + ".csv";

        try(BufferedWriter writer = new BufferedWriter(new OutputStreamWriter(
                new FileOutputStream(save_path), StandardCharsets.UTF_8))) {
            Molecule mol;
            RxnMolecule reaction;
            RxnMolecule reactionCore;
            RxnMolecule reaction1Neighbor;
            Standardizer std = new Standardizer(new File("./strip_salts.xml"));
            for (int i = 0; i < molStreams.size(); i++) {
                mol = getMolFromInputStream(molStreams.get(i));
                reaction = RxnMolecule.getReaction(mol);
                String sMol = MolExporter.exportToFormat(reaction, FORMAT);
                reaction = RxnMolecule.getReaction(MolImporter.importMol(sMol));
                stripSalts(reaction, std);
                if (checkMwOfProductsInReaction(reaction)) {
                    continue;
                }
                if (reaction.getProductCount() != 1 | reaction.getReactantCount() > 3) {
                    continue;
                }
                Molecule product = reaction.getProduct(0);
                mapReaction(reaction, MappingStyle.CHANGING);
                RxnMolecule reactionCloneForCore = reaction.clone();
                RxnMolecule reactionCloneFor1Neighbor = reaction.clone();
                getReactionTemplate(reactionCloneForCore, "0");
                getReactionTemplate(reactionCloneFor1Neighbor, "1");
                if (reactionCloneForCore.getProduct(0).isEmpty() | reactionCloneFor1Neighbor.getProduct(0).isEmpty()) {
                    continue;
                }
                AutoMapper.unmap(product);
                AutoMapper.unmap(reaction);
                AutoMapper.unmap(reactionCloneForCore);
                AutoMapper.unmap(reactionCloneFor1Neighbor);
                reactionCloneForCore.removeEmptyComponents();
                reactionCloneFor1Neighbor.removeEmptyComponents();

                String rCore = MolExporter.exportToFormat(reactionCloneForCore, FORMAT);
                reactionCore = RxnMolecule.getReaction(MolImporter.importMol(rCore));
                String r1Neighbor = MolExporter.exportToFormat(reactionCloneFor1Neighbor, FORMAT);
                reaction1Neighbor = RxnMolecule.getReaction(MolImporter.importMol(r1Neighbor));
                if (reactionCore.getProductCount() != 1 | reaction1Neighbor.getProductCount() != 1) {
                    continue;
                }
                writer.write(idList.get(i) + "," +
                        MolExporter.exportToFormat(product, FORMAT) + "," +
                        MolExporter.exportToFormat(sortReactantsInReaction(reaction), FORMAT) + "," +
                        MolExporter.exportToFormat(sortReactantsInReaction(reactionCore), FORMAT) + "," +
                        MolExporter.exportToFormat(sortReactantsInReaction(reaction1Neighbor), FORMAT) + "," +
                        max_pub_year);
                writer.newLine();
            }
            System.out.println("[SAVE] " + save_path);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}

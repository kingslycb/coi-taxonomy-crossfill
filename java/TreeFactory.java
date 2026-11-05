/*
 * TreeFactory.java
 *
 * Copyright 2006 Michigan State University Board of Trustees
 * Created on June 24, 2002, 5:01 PM
 *
 * From the RDP Classifier 2.14 (https://github.com/rdpstaff/classifier), distributed under the GNU
 * General Public License, version 2 or (at your option) any later version. Modified for this workflow:
 *  1. The word-prior and word-pointer arrays are sized from the word length set in GoodWordIterator
 *     (the original fixes them at 65,536 words, i.e. word length 8); needed for word length 7.
 *  2. The duplicate-name check counts (name, rank, parent) instead of (name, rank).
 * Compiled ahead of classifier.jar on the class path, it replaces the original class.
 */
package edu.msu.cme.rdp.classifier.train;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Reader;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.StringTokenizer;

/**
 * A TreeFactory reads the raw taxonomic information and creates the RawHierarchyTree nodes.
 * @author  wangqion
 * @version MODIFIED: homonym (name, rank) duplicates are allowed when they occur under
 * DIFFERENT parents. The classifier already disambiguates by (name, parentID, depth)
 * inside getTaxonomy(); the original global sanity check rejected legitimate NCBI homonyms.
 */
public class TreeFactory {

    private static String dformat = "%1$.2f";
    private List genusNodeList = null;
    private List genus_wordConditionalProbList = new ArrayList();
    private int[] wordProbPointerArr;
    private float[] logArr;
    private RawHierarchyTree rootTree;
    private Map<String, List> taxnameMap = new HashMap<String, List>();
    private Map<Integer, Taxonomy> taxidMap = new HashMap<Integer, Taxonomy>();
    private Map<String, Integer> taxnameRankMap = new HashMap<String, Integer>();
    private float[] wordPriorArr;
    private int ROOT_DEPTH = 0;
    private int totalSequences = 0;
    private final float WF1 = (float) 0.5;
    private final float WF2 = (float) 1;
    private BufferedWriter treeFile;
    private String trainingVersion;
    private String trained_rank = null;

    public TreeFactory(Reader taxReader, int trainsetNo, String version, String modification) throws IOException, NameRankDupException {
        edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator gw =
            new edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator("ACGT");
        int numWords = gw.getMask() + 1;
        wordPriorArr = new float[numWords];
        wordProbPointerArr = new int[numWords + 1];
        trainingVersion = "<trainsetNo>" + trainsetNo + "</trainsetNo>" + "<version>" + version + "</version><modversion>" + modification + "</modversion>";
        creatTaxidMap(taxReader);
    }

    /** It reads in a file containing the taxonomy information for all the nodes.
     * The taxonomy format is taxid(int), taxname(string), parentid(int),
     * depth(int) and hierarchy level(string) seperated by * in one line.
     */
    private void creatTaxidMap(Reader taxReader) throws IOException, NameRankDupException {
        BufferedReader reader = new BufferedReader(taxReader);
        String line;

        while ((line = reader.readLine()) != null) {
            if (line.length() == 0) {
                continue;
            }
            StringTokenizer st = new StringTokenizer(line, "*");
            if (st.countTokens() < 5) {
                throw new IllegalArgumentException("\nIllegal taxonomy format at " + line);
            }
            try {
                int taxid = Integer.parseInt(st.nextToken().trim());
                String taxname = st.nextToken().trim();
                int pid = Integer.parseInt(st.nextToken().trim());
                int depth = Integer.parseInt(st.nextToken().trim());
                List taxList = (ArrayList) taxnameMap.get(taxname);
                if (taxList == null) {
                    taxList = new ArrayList();
                }
                Taxonomy tax = new Taxonomy(taxid, taxname, pid, depth, st.nextToken().trim());
                taxList.add(tax);
                taxnameMap.put(taxname, taxList);
                if ( taxidMap.containsKey(taxid)){
                    throw new NameRankDupException("Error: duplicate taxid found : " + taxid);
                }
                taxidMap.put(new Integer(taxid), tax);
                String name_rank_parent = (taxname + "\t" + tax.hierLevel + "\t" + pid).toLowerCase();
                Integer existCount = taxnameRankMap.get(name_rank_parent);
                if (existCount == null) {
                    existCount = new Integer(0);
                }
                taxnameRankMap.put(name_rank_parent, existCount.intValue() + 1);
                // check the root tree
                if (tax.depth == this.ROOT_DEPTH) {
                    if (rootTree == null) {
                        rootTree = new RawHierarchyTree(taxname, null, tax);
                    } else {
                        throw new IllegalArgumentException("Error: taxon " + tax.getTaxID() + " has the depth set to '0'. Only the root taxon can have the depth set to '0'");
                    }
                }
            } catch (NumberFormatException e) {
                throw new IllegalArgumentException("\nError: "
                        + "The value for taxid, parentid and the depth should be integer in : " + line);
            }
        }

        if (rootTree == null) {
            throw new IllegalArgumentException("Error: no root taxon with depth '0' defined in the taxonomy file.");
        }
        String errors = nameRankSanityCheck();
        if (errors != null) {
            throw new NameRankDupException("Error: duplicate taxon name and rank under the same parent in the taxonomy file.\n" + errors);
        }
    }

    private String nameRankSanityCheck() {
        StringBuilder buf = new StringBuilder();
        for (String name_rank_parent : taxnameRankMap.keySet()) {
            int count = taxnameRankMap.get(name_rank_parent);
            if (count > 1) {
                buf.append(name_rank_parent).append("\t").append(count).append("\n");
            }
        }
        if (buf.length() == 0) {
            return null;
        }
        return buf.toString();
    }

    public void parseSequenceFile(LineageSequenceParser parser) throws IOException {
        while (parser.hasNext()) {
            LineageSequence pSeq = parser.next();
            if (pSeq.getAncestors().size() == 1) {  // only contains the taxid of the immediate taxon
                addSequencewithTaxid(pSeq);
            } else {
                addSequencewithLineage(pSeq);
            }
        }
        parser.close();
    }

    public void parseCopyNumberFile(String cnFile) throws IOException {
         HashMap<String, CopyNumberParser.TaxonCopyNumber> cnMap = new CopyNumberParser().parse(cnFile);
         this.getGenusNodeList();
         for (int index = 0; index < genusNodeList.size(); index++) {
                RawHierarchyTree aTree = ((RawHierarchyTree) genusNodeList.get(index));
                String key = (aTree.getName() + aTree.getTaxonomy().getHierLevel()).toLowerCase();
                CopyNumberParser.TaxonCopyNumber cn = cnMap.get(key);
                if ( cn != null){
                    aTree.setCopyNumber(cn.getCopyNumber());
                    cnMap.remove(key);
                }
         }
         StringBuilder msg = new StringBuilder();
         for ( CopyNumberParser.TaxonCopyNumber cn: cnMap.values()){
             if ( !cn.rank.equalsIgnoreCase(this.trained_rank)) continue;
             msg.append(cn.name).append("\n");
         }
         if (msg.length() > 0){
             throw new IllegalArgumentException("Error: The following taxon names at rank " + this.trained_rank + " in the copynumber file are not found in the input taxonomy file\n" + msg);
         }
         bottomUpSetCopyNumber(this.rootTree);
         if ( this.rootTree.hasCopyNumber()){
            topDownFillCopyNumber(this.rootTree);
         } else {
             throw new IllegalArgumentException("Error: Problem setting copy number, no taxon at rank " + this.trained_rank + " found in the copynumber file");
         }
    }

    private void addSequencewithTaxid(LineageSequence pSeq) throws IOException {
        ArrayList<Taxonomy> taxonList = new ArrayList<Taxonomy>();
        Taxonomy tax = this.taxidMap.get(new Integer((String) pSeq.getAncestors().get(0)));

        while (tax != null) {
            taxonList.add(tax);
            if (tax.getTaxID() == rootTree.getTaxonomy().getTaxID()) {
                break;
            }
            tax = this.taxidMap.get(new Integer(tax.getParentID()));
        }

        if (tax == null) {
            throw new IllegalArgumentException("Problem retrieving ancestor taxon for Sequence " + pSeq.getSeqName());
        }
        if (tax.getTaxID() != rootTree.getTaxonomy().getTaxID()) {
            throw new IllegalArgumentException("Sequence " + pSeq.getSeqName()
                    + " has conflicting ancestor root name: " + tax.getTaxName());
        }
        RawHierarchyTree curTree = rootTree;
        for (int i = taxonList.size() - 2; i >= 0; i--) {
            RawHierarchyTree tmp = curTree.getSubclassbyName(taxonList.get(i).getTaxName());
            if (tmp == null) {
                curTree = new RawHierarchyTree(taxonList.get(i).getTaxName(), curTree, taxonList.get(i));
            } else {
                curTree = tmp;
            }
            if (i == 0) {
                curTree.initWordOccurrence(pSeq, wordPriorArr);
                if (this.trained_rank == null) {
                    this.trained_rank = curTree.getTaxonomy().getHierLevel();
                } else if (!this.trained_rank.equalsIgnoreCase(curTree.getTaxonomy().getHierLevel())) {
                    throw new IllegalArgumentException("Sequence " + pSeq.getSeqName()
                            + " has different lowest rank: " + curTree.getTaxonomy().getHierLevel() + " from the previous lowest rank: " + this.trained_rank);
                }
                totalSequences++;
            }
        }
    }

    private void addSequencewithLineage(LineageSequence pSeq) throws IOException {
        int size = pSeq.getAncestors().size();
        if (!((String) pSeq.getAncestors().get(0)).equalsIgnoreCase(rootTree.getName())) {
            throw new IllegalArgumentException("Sequence " + pSeq.getSeqName()
                    + " has conflicting root name: " + pSeq.getAncestors().get(0));
        }
        RawHierarchyTree curTree = rootTree;
        for (int i = 1; i < size; i++) {
            RawHierarchyTree tmp = curTree.getSubclassbyName(pSeq.getAncestors().get(i));
            if (tmp == null) {
                Taxonomy tax = getTaxonomy(pSeq, ((Taxonomy) curTree.getTaxonomy()).taxID, i);
                curTree = new RawHierarchyTree((String) pSeq.getAncestors().get(i), curTree, tax);
            } else {
                curTree = tmp;
            }
            if (i == size - 1) {
                curTree.initWordOccurrence(pSeq, wordPriorArr);
                totalSequences++;
                if (this.trained_rank == null) {
                    this.trained_rank = curTree.getTaxonomy().getHierLevel();
                } else if (!this.trained_rank.equalsIgnoreCase(curTree.getTaxonomy().getHierLevel())) {
                    throw new IllegalArgumentException("Sequence " + pSeq.getSeqName()
                            + " has different lowest rank: " + curTree.getTaxonomy().getHierLevel() + " from the previous lowest rank: " + this.trained_rank);
                }
            }
        }
    }

    private Taxonomy getTaxonomy(LineageSequence pSeq, int pid, int index) {
        List ancestor = pSeq.getAncestors();
        if (ancestor.isEmpty()) {
            throw new IllegalArgumentException("Error: No ancestors found for sequence: " + pSeq.getSeqName()
                    + "! Please check the source file.");
        }
        String name = (String) ancestor.get(index);
        List taxList = (ArrayList) taxnameMap.get(name);
        if (taxList == null) {
            throw new IllegalArgumentException("\nThe taxID for ancestor: " + name + " of sequence: " + pSeq.getSeqName()
                    + " at depth: " + index + " with parent id: " + pid + " is not found!");
        }

        Taxonomy result = null;
        for (int i = 0; i < taxList.size(); i++) {
            Taxonomy tax = (Taxonomy) taxList.get(i);
            if (tax.parentID == pid && tax.depth == index) {
                result = tax;
                break;
            }
        }
        if (result == null) {
            throw new IllegalArgumentException("\nThe taxID for ancestor: " + name + " of sequence: " + pSeq.getSeqName()
                    + " at depth: " + index + " with parent id: " + pid + " is not found!");
        }
        return result;
    }

    public RawHierarchyTree getRoot() {
        return rootTree;
    }

    void createGenusWordConditionalProb() {
        getGenusNodeList();
        if (genusNodeList.isEmpty()) {
            throw new IllegalArgumentException("\nThere is no node at " + this.trained_rank);
        }
        int maxNumOfLeaves = 0;

        for (int i = 0; i < wordPriorArr.length; i++) {
            wordPriorArr[i] = (wordPriorArr[i] + WF1) / (totalSequences + WF2);
            wordProbPointerArr[i] = genus_wordConditionalProbList.size();
            for (int index = 0; index < genusNodeList.size(); index++) {
                RawHierarchyTree aTree = ((RawHierarchyTree) genusNodeList.get(index));
                int wordOccurrence = aTree.getWordOccurrence(i);
                int numOfLeaves = aTree.getLeaveCount();
                if (wordOccurrence > 0) {
                    float prob = (float) Math.log((wordOccurrence + wordPriorArr[i]) / (numOfLeaves + WF2));
                    genus_wordConditionalProbList.add(new RawGenusWordConditionalProb(index, prob));
                }
                if (numOfLeaves > maxNumOfLeaves) {
                    maxNumOfLeaves = numOfLeaves;
                }
            }
            wordPriorArr[i] = (float) Math.log(wordPriorArr[i]);
        }

        wordProbPointerArr[wordProbPointerArr.length - 1] = genus_wordConditionalProbList.size();

        for (int index = 0; index < genusNodeList.size(); index++) {
            ((RawHierarchyTree) genusNodeList.get(index)).releaseWordOccurrence();
        }

        logArr = new float[++maxNumOfLeaves];
        for (int i = 0; i < maxNumOfLeaves; i++) {
            logArr[i] = (float) Math.log(i + WF2);
        }
    }

    float getLogWordPrior(int wordIndex) {
        return wordPriorArr[wordIndex];
    }

    List getGenusNodeList() {
        if (this.trained_rank == null) {
            throw new IllegalArgumentException("Need to read a sequence file to set the rank to train on");
        }
        if ( genusNodeList != null) return genusNodeList;
        genusNodeList = new ArrayList();
        createNodeList(this.getRoot(), this.trained_rank, genusNodeList);
        return genusNodeList;
    }

    float getLogLeaveCount(int i) {
        return logArr[i];
    }

    int getStartIndex(int wordIndex) {
        return wordProbPointerArr[wordIndex];
    }

    int getStopIndex(int wordIndex) {
        return wordProbPointerArr[wordIndex + 1];
    }

    RawGenusWordConditionalProb getWordConditionalProb(int posIndex) {
        return (RawGenusWordConditionalProb) genus_wordConditionalProbList.get(posIndex);
    }

    void createNodeList(RawHierarchyTree node, String level, List nodeList) {
        if (node == null) {
            return;
        }
        if (((Taxonomy) node.getTaxonomy()).hierLevel.equalsIgnoreCase(level)) {
            nodeList.add(node);
            node.setGenusIndex(nodeList.size() - 1);
            return;
        }
        Collection al = new ArrayList();
        if ((al = node.getSubclasses()).isEmpty()) {
            return;
        }
        Iterator i = al.iterator();
        while (i.hasNext()) {
            createNodeList((RawHierarchyTree) i.next(), level, nodeList);
        }
    }

    void bottomUpSetCopyNumber(RawHierarchyTree node) {
        if ( node.getTaxonomy().hierLevel.equals(this.trained_rank)) return;
        Collection<RawHierarchyTree> al = node.getSubclasses();
        if (al.isEmpty())  return;
        Iterator i = al.iterator();
        while (i.hasNext()) {
            bottomUpSetCopyNumber((RawHierarchyTree) i.next());
        }
        float sum = 0f;
        int childwithcn = 0;
        for (RawHierarchyTree c: al){
            if ( c.hasCopyNumber()){
                sum += c.getCopyNumber();
                childwithcn ++;
            }
        }
        node.setCopyNumber(sum/childwithcn);
    }

    void topDownFillCopyNumber(RawHierarchyTree node) {
        if ( !node.hasCopyNumber()){
            node.setCopyNumber(node.getParent().getCopyNumber());
        }
        Collection<RawHierarchyTree> al = node.getSubclasses();
        if (al.isEmpty())     return;
        Iterator i = al.iterator();
        while (i.hasNext()) {
            topDownFillCopyNumber((RawHierarchyTree) i.next());
        }
    }

    void printTrainingFiles(String outdir) throws IOException {
        treeFile = new BufferedWriter(new FileWriter(outdir + "bergeyTrainingTree.xml"));
        treeFile.write(trainingVersion + "<file>bergeyTrainingTree</file>\n");
        displayTrainingTree(this.rootTree);
        treeFile.close();
    }

    private void displayTrainingTree(RawHierarchyTree root) throws IOException {
        Taxonomy taxon = ((Taxonomy) root.getTaxonomy());
        treeFile.write("<TreeNode name=\"" + root.getName().replaceAll("&", "").replaceAll("\"", "&quot;") + "\" taxid=\""
                + taxon.taxID + "\" rank=\"" + taxon.hierLevel + "\" parentTaxid=\""
                + taxon.parentID + "\" leaveCount=\""
                + root.getLeaveCount() + "\" genusIndex=\"" + root.getGenusIndex() + "\"");
        if (rootTree.hasCopyNumber()){
            treeFile.write(" cpNumber=\"" + String.format(dformat, root.getCopyNumber()) + "\"");
        }
        treeFile.write("></TreeNode>\n");
        Iterator i = root.getSubclasses().iterator();
        while (i.hasNext()) {
            displayTrainingTree((RawHierarchyTree) i.next());
        }
    }

    void printWordPriors(String outdir) throws IOException {
        BufferedWriter outfile = new BufferedWriter(new FileWriter(outdir + "logWordPrior.txt"));
        outfile.write(trainingVersion + "<file>logWordPrior</file>\n");
        for (int i = 0; i < wordPriorArr.length; i++) {
            outfile.write(i + "\t" + wordPriorArr[i] + "\n");
        }
        outfile.close();
    }

    void printWordConditionalProbIndexArr(String outdir) throws IOException {
        BufferedWriter outfile = new BufferedWriter(new FileWriter(outdir + "wordConditionalProbIndexArr.txt"));
        outfile.write(trainingVersion + "<file>wordConditionalProbIndexArr</file>\n");
        for (int i = 0; i < wordProbPointerArr.length; i++) {
            outfile.write(i + "\t" + wordProbPointerArr[i] + "\n");
        }
        outfile.close();
    }

    void printGenusIndex_WordProbArr(String outdir) throws IOException {
        BufferedWriter outfile = new BufferedWriter(new FileWriter(outdir + "genus_wordConditionalProbList.txt"));
        outfile.write(trainingVersion + "<file>genus_wordConditionalProbList</file>\n");
        for (int i = 0; i < genus_wordConditionalProbList.size(); i++) {
            RawGenusWordConditionalProb prob = (RawGenusWordConditionalProb) genus_wordConditionalProbList.get(i);
            outfile.write(prob.getGenusIndex() + "\t" + prob.getProbability() + "\n");
        }
        outfile.close();
    }
}
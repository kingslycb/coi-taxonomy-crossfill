/*
 * WordSizeMain - runs the RDP Classifier (train or classify) with a word length other than the default 8.
 * The word length is read from the system property 'wordsize' (default 7), e.g.
 *   java -Dwordsize=7 -cp classes:classifier.jar:lib/* WordSizeMain train -o model -s seqs.fasta -t taxonomy.txt
 *   java -Dwordsize=7 -cp classes:classifier.jar:lib/* WordSizeMain classify -t model/rRNAClassifier.properties ...
 * The same word length must be used for training and classification.
 */
import edu.msu.cme.rdp.classifier.cli.ClassifierMain;
import edu.msu.cme.rdp.readseq.utils.orientation.GoodWordIterator;

public class WordSizeMain {
    public static void main(String[] args) throws Exception {
        GoodWordIterator.setWordSize(Integer.parseInt(System.getProperty("wordsize", "7")));
        ClassifierMain.main(args);
    }
}

/*
 * MIT License
 *
 * Copyright (c) 2020 Jonas Schaub, Achim Zielesny, Christoph Steinbeck, Maria Sorokina
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package de.unijena.cheminf.deglycosylation.stats;

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
import com.mongodb.client.MongoClient;
import com.mongodb.client.MongoClients;
import com.mongodb.client.MongoCollection;
import com.mongodb.client.MongoCursor;
import com.mongodb.client.MongoDatabase;
import de.unijena.cheminf.deglycosylation.SugarRemovalUtility;
import org.bson.Document;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.awt.Color;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.Arrays;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Locale;
import java.util.Objects;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * This test class uses the <a href="https://doi.org/10.1186/s13321-020-00467-y"
 * >"Sugar Removal Utility (SRU)"</a> to detect glycosylation in various databases, mainly the <a href="https://doi.org/10.1186/s13321-020-00478-9"
 * >"COlleCtion of Open Natural prodUcTs (COCONUT)"</a>. It is imported via MongoDB. Therefore, MongoDB needs to be running
 * on the platform. See credentials for making the connection in the private static final constants.
 * Also, most of the resource files (different data sets) are not included in the repository.
 * Results obtained with these analyses are published in <a href="https://doi.org/10.3390/biom11040486"
 * >"Description and Analysis of Glycosidic Residues in the Largest Open Natural Products Database (Schaub et al. 2021)"</a>
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class GlycosylationStatisticsTest extends SugarRemovalUtility {
    //<editor-fold desc="Private static final constants">
    /**
     * Host of the MongoDB instance
     */
    private static final String HOST = "localhost";

    /**
     * Port where the MongoDB instance is running
     */
    private static final int PORT = 27017;

    /**
     * Name of the MongoDB database to access; this COCONUT version contains 401,624 unique natural products
     */
    private static final String DATABASE_NAME = "COCONUT2020november03";

    /**
     * Collection from the database to load
     */
    private static final String COLLECTION_NAME = "uniqueNaturalProduct";

    /**
     * Name of the COCONUT SDF to use
     */
    private static final String SDF_NAME = "COCONUT_DB_november_18.sdf";

    /**
     * Name of the output folder for this test class
     */
    private static final String OUTPUT_FOLDER_NAME = "GlycosylationStatisticsTest_Output";

    /**
     * Name of the document variable that contains an ID of the given molecule, "coconut_id" or "inchikey"
     */
    private static final String ID_KEY = "coconut_id";

    /**
     * Name of the document variable that contains the SMILES code string of the given molecule, "clean_smiles" or
     * "smiles" or "unique_smiles"
     */
    private static final String SMILES_CODE_KEY = "unique_smiles";

    /**
     * Name of the ZINC 'biogenic' subset SMILES file in the resource directory
     */
    private static final String ZINC_BIOGENIC_SUBSET_FILE_NAME = "ZINC_biogenic_subset_2020_Okt_17.txt";

    /**
     * Separator for created CSV output files
     */
    private static final String OUTPUT_FILE_SEPARATOR = ";";

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(GlycosylationStatisticsTest.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Test methods">
    //<editor-fold desc="Tests involving databases">
    //<editor-fold desc="COCONUT">
    /**
     * This test method visualises deglycosylation results of all sugar-containing molecules in COCONUT. For each of these
     * molecules, an image file depicting the original structure and one depicting the deglycosylated molecule (aglycon)
     * are created in the respective output directory (./GlycosylationStatisticsTest_Output/coconut_deglycosylation_visualization_test/).
     * To run this test, remove the JUnit @Ignore tag that is set because
     * this test creates many files. Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void coconutDeglycosylationVisualizationTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_deglycosylation_visualization_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //for depiction without explicit hydrogens:
                tmpMolecule = AtomContainerManipulator.removeHydrogens(tmpMolecule);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not detected/removed/counted!
                // note also: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                boolean tmpHasAnyTypeOfSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpHasAnyTypeOfSugarsCounter++;
                    tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath
                            + File.separator + tmpID + "_deglycosylated.png");
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println();
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        System.out.println(tmpHasAnyTypeOfSugarsCounter * 2 + " image files were created.");
        tmpCursor.close();
    }

    /**
     * This test method does some basic glycosylation statistics of COCONUT. It calculates: How many molecules have sugars,
     * circular sugars, linear sugars, both, and how many molecules are basically sugars. All statistics are printed to
     * console and also compiled in an output file created in the directory ./GlycosylationStatisticsTest_Output/coconut_stats_basics_test/.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsBasicsTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_basics_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsCNPs = new ArrayList<>(50000);
        int tmpHasNoSugarsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList<>(2000);
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList<>(700);
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarCNPs = new ArrayList<>(2000);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not detected/removed/counted!
                // note also: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                //note: this property will also be true if sugars were detected but not removed because they are not terminal
                boolean tmpHasAnyTypeOfSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyCircularSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyLinearSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpHasAnyTypeOfSugarsCounter++;
                    tmpHasAnyTypeOfSugarsCNPs.add(tmpID);
                    if (tmpHasAnyCircularSugar) {
                        tmpHasCircularSugarsCounter++;
                        tmpHasCircularSugarsCNPs.add(tmpID);
                    }
                    if (tmpHasAnyLinearSugar) {
                        tmpHasLinearSugarsCounter++;
                        tmpHasLinearSugarsCNPs.add(tmpID);
                    }
                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsCNPs.add(tmpID);
                    }
                    if (tmpDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                        tmpBasicallyASugarCNPs.add(tmpID);
                    }
                } else {
                    tmpHasNoSugarsCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        tmpOutputWriter.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        double tmpPercentage = ((double) tmpHasAnyTypeOfSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain sugars.");
        System.out.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        tmpOutputWriter.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpOutputWriter.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        tmpOutputWriter.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Sugar-containing molecules: " + tmpHasAnyTypeOfSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Circular-sugar-containing molecules: " + tmpHasCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Linear-sugar-containing molecules: " + tmpHasLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules containing both circular and linear sugars: " + tmpHasCircularAndLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a sugar: " + tmpBasicallyASugarCNPs);
        tmpOutputWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpMoleculesCounter, tmpHasNoSugarsCounter + tmpHasAnyTypeOfSugarsCounter);
    }

    /**
     * This test method does some glycosylation statistics focused on circular sugars of COCONUT. It calculates:
     * How many molecules have circular sugars, only circular sugars (no linear ones), terminal circular sugars,
     * only terminal circular sugars (no non-terminal circular ones), non-terminal circular sugars, only non-terminal
     * circular sugars (no terminal circular ones), terminal and non-terminal circular sugars, circular sugars with
     * glycosidic bonds, terminal circular sugars with glycosidic bonds, non-terminal circular sugars with glycosidic
     * bonds, both, spiro sugars, and how many molecules qualify for the SRU glycosidic bond exemption. All statistics
     * are printed to console and also compiled in an output file created in the directory ./GlycosylationStatisticsTest_Output/coconut_stats_circular_sugars_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsCircularSugarsTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_circular_sugars_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalCircularSugarsCNPs = new ArrayList<>(40000);
        int tmpHasNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasNonTerminalCircularSugarsCNPs = new ArrayList<>(12000);
        int tmpHasOnlyTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalCircularSugarsCNPs = new ArrayList<>(40000);
        int tmpHasOnlyNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalCircularSugarsCNPs = new ArrayList<>(10000);
        int tmpHasTerminalAndNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalCircularSugarsCNPs = new ArrayList<>(5000);
        int tmpHasOnlyCircularSugarsCounter = 0;
        List<String> tmpHasOnlyCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasGlycosidicBondCounter = 0;
        List<String> tmpHasGlycosidicBondCNPs = new ArrayList<>(50000);
        int tmpHasGlycosidicBondOnTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnTerminalSugarCNPs = new ArrayList<>(40000);
        int tmpHasGlycosidicBondOnNonTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnNonTerminalSugarCNPs = new ArrayList<>(12000);
        int tmpGlycosidicBondExemptionCounter = 0;
        List<String> tmpGlycosidicBondExemptionCNPs = new ArrayList<>(100);
        int tmpHasSpiroSugarCounter = 0;
        List<String> tmpHasSpiroSugarCNPs = new ArrayList<>(1000);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpSugarRemovalUtil.hasCircularSugars(tmpMolecule);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpSugarRemovalUtil.hasLinearSugars(tmpMolecule);
                if (tmpHasAnyCircularSugar) {
                    tmpHasCircularSugarsCounter++;
                    tmpHasCircularSugarsCNPs.add(tmpID);
                    if (!tmpHasAnyLinearSugar) {
                        tmpHasOnlyCircularSugarsCounter++;
                        tmpHasOnlyCircularSugarsCNPs.add(tmpID);
                    }
                    //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                    List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                    int tmpNumberOfCircularSugarMoieties;
                    int tmpNumberOfTerminalCircularSugarMoieties;
                    int tmpNumberOfNonTerminalCircularSugarMoieties;
                    int tmpNumberOfGlycosidicBonds;
                    int tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                    int tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond;
                    tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();
                    //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                    List<IAtomContainer> tmpRemovedTerminalCircularSugarMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                    //-1 for the deglycosylated core at the beginning of the list
                    tmpNumberOfTerminalCircularSugarMoieties = tmpRemovedTerminalCircularSugarMoieties.size() - 1 ;
                    tmpNumberOfNonTerminalCircularSugarMoieties = tmpNumberOfCircularSugarMoieties - tmpNumberOfTerminalCircularSugarMoieties;
                    Assert.assertTrue(tmpNumberOfNonTerminalCircularSugarMoieties >= 0);
                    //leaving default! Now, only circular sugars having glycosidic bonds are in the candidates and removed moieties
                    tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
                    boolean tmpMoleculeQualifiesForExemption = tmpSugarRemovalUtil.isQualifiedForGlycosidicBondExemption(tmpMolecule.clone());
                    if (tmpMoleculeQualifiesForExemption) {
                        //note: these molecules are basically made up of one circular sugar without a glycosidic bond,
                        // so that the whole molecule is empty after the removal of this sugar moiety. This is an
                        // exemption implemented in the compilation of circular sugar moieties with detection of glycosidic bonds.
                        tmpGlycosidicBondExemptionCounter++;
                        tmpGlycosidicBondExemptionCNPs.add(tmpID);
                        tmpNumberOfGlycosidicBonds = 0;
                        tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = 0;
                        tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond = 0;
                    } else {
                        List<IAtomContainer> tmpCircularSugarCandidatesWithGlycosidicBondsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                        tmpNumberOfGlycosidicBonds = tmpCircularSugarCandidatesWithGlycosidicBondsList.size();
                        List<IAtomContainer> tmpRemovedTerminalCircularMoietiesWithGlycosidicBond = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = tmpRemovedTerminalCircularMoietiesWithGlycosidicBond.size() - 1;
                        tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond = tmpNumberOfGlycosidicBonds - tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                        Assert.assertTrue(tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond >= 0);
                    }
                    //back to default!
                    tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
                    if (tmpNumberOfTerminalCircularSugarMoieties > 0) {
                        tmpHasTerminalCircularSugarsCounter++;
                        tmpHasTerminalCircularSugarsCNPs.add(tmpID);
                        if (tmpNumberOfNonTerminalCircularSugarMoieties == 0) {
                            tmpHasOnlyTerminalCircularSugarsCounter++;
                            tmpHasOnlyTerminalCircularSugarsCNPs.add(tmpID);
                        }
                    }
                    if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                        tmpHasNonTerminalCircularSugarsCounter++;
                        tmpHasNonTerminalCircularSugarsCNPs.add(tmpID);
                        if (tmpNumberOfTerminalCircularSugarMoieties == 0) {
                            tmpHasOnlyNonTerminalCircularSugarsCounter++;
                            tmpHasOnlyNonTerminalCircularSugarsCNPs.add(tmpID);
                        }
                    }
                    if (tmpNumberOfTerminalCircularSugarMoieties > 0 && tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                        tmpHasTerminalAndNonTerminalCircularSugarsCounter++;
                        tmpHasTerminalAndNonTerminalCircularSugarsCNPs.add(tmpID);
                    }
                    if (tmpNumberOfGlycosidicBonds > 0) {
                        tmpHasGlycosidicBondCounter++;
                        tmpHasGlycosidicBondCNPs.add(tmpID);
                        if (tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond > 0) {
                            tmpHasGlycosidicBondOnTerminalSugarCounter++;
                            tmpHasGlycosidicBondOnTerminalSugarCNPs.add(tmpID);
                        }
                        if (tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond > 0) {
                            tmpHasGlycosidicBondOnNonTerminalSugarCounter++;
                            tmpHasGlycosidicBondOnNonTerminalSugarCNPs.add(tmpID);
                        }
                    }
                }
                tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(true);
                boolean tmpHasCircularSugarsIncludingSpiro = tmpSugarRemovalUtil.hasCircularSugars(tmpMolecule);
                if (tmpHasCircularSugarsIncludingSpiro) {
                    int tmpCircularCandidateNumberInlcudingSpiro = tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpMolecule);
                    tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);
                    int tmpCircularCandidateNumberExcludingSpiro = tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpMolecule);
                    Assert.assertTrue(tmpCircularCandidateNumberExcludingSpiro <= tmpCircularCandidateNumberInlcudingSpiro);
                    if (tmpCircularCandidateNumberInlcudingSpiro > tmpCircularCandidateNumberExcludingSpiro) {
                        tmpHasSpiroSugarCounter++;
                        tmpHasSpiroSugarCNPs.add(tmpID);
                    }
                }
                tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(false);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpOutputWriter.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        double tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Only circular sugar containing molecules counter: " + tmpHasOnlyCircularSugarsCounter);
        tmpOutputWriter.println("Only circular sugar containing molecules counter: " + tmpHasOnlyCircularSugarsCounter);
        System.out.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasOnlyCircularSugarsCounter) + " molecules " +
                " contain circular sugars and also linear sugar moieties.");
        tmpOutputWriter.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasOnlyCircularSugarsCounter) + " molecules " +
                "contain circular sugars and also linear sugar moieties.");
        System.out.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        System.out.println("Only terminal circular sugar containing molecules counter: " + tmpHasOnlyTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Only terminal circular sugar containing molecules counter: " + tmpHasOnlyTerminalCircularSugarsCounter);
        System.out.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        System.out.println("Only non-terminal circular sugar containing molecules counter: " + tmpHasOnlyNonTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Only non-terminal circular sugar containing molecules counter: " + tmpHasOnlyNonTerminalCircularSugarsCounter);
        System.out.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        System.out.println("Circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondCounter);
        tmpOutputWriter.println("Circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondCounter);
        System.out.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasGlycosidicBondCounter) + " molecules " +
                "only have circular sugar moieties that are NOT attached via a glycosidic bond.");
        tmpOutputWriter.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasGlycosidicBondCounter) + " molecules " +
                "only have circular sugar moieties that are NOT attached via a glycosidic bond.");
        System.out.println("Terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnTerminalSugarCounter);
        tmpOutputWriter.println("Terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on non-terminal circular sugar moieties.");
        tmpOutputWriter.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on non-terminal circular sugar moieties.");
        System.out.println("Non-terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnNonTerminalSugarCounter);
        tmpOutputWriter.println("Non-terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnNonTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnNonTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on terminal circular sugar moieties.");
        tmpOutputWriter.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnNonTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on terminal circular sugar moieties.");
        int tmpHasBoth = tmpHasGlycosidicBondOnNonTerminalSugarCounter + tmpHasGlycosidicBondOnTerminalSugarCounter - tmpHasGlycosidicBondCounter;
        System.out.println(tmpHasBoth + " molecules have both, terminal and non-terminal sugar moieties attached via a glycosidic bond.");
        tmpOutputWriter.println(tmpHasBoth + " molecules have both, terminal and non-terminal sugar moieties attached via a glycosidic bond.");
        System.out.println("Molecules that qualify for the glycosidic bond exemption counter: " + tmpGlycosidicBondExemptionCounter);
        tmpOutputWriter.println("Molecules that qualify for the glycosidic bond exemption counter: " + tmpGlycosidicBondExemptionCounter);
        System.out.println("Molecules that have spiro sugars counter: " + tmpHasSpiroSugarCounter);
        tmpOutputWriter.println("Molecules that have spiro sugars counter: " + tmpHasSpiroSugarCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Circular sugar containing molecules: " + tmpHasCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal circular sugar containing molecules: " + tmpHasTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Non-terminal circular sugar containing molecules: " + tmpHasNonTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only terminal circular sugar containing molecules: " + tmpHasOnlyTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only non-terminal circular sugar containing molecules: " + tmpHasOnlyNonTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal and non-terminal circular sugar containing molecules: " + tmpHasTerminalAndNonTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only circular sugar containing molecules: " + tmpHasOnlyCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Circular sugar with glycosidic bond containing molecules: " + tmpHasGlycosidicBondCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal circular sugar with glycosidic bond containing molecules: " + tmpHasGlycosidicBondOnTerminalSugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Non-terminal circular sugar with glycosidic bond containing molecules: " + tmpHasGlycosidicBondOnNonTerminalSugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules that qualify for the glycosidic bond exemption: " + tmpGlycosidicBondExemptionCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules that have spiro sugars: " + tmpHasSpiroSugarCNPs);
        tmpOutputWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpHasTerminalAndNonTerminalCircularSugarsCounter
                + tmpHasOnlyTerminalCircularSugarsCounter
                + tmpHasOnlyNonTerminalCircularSugarsCounter);
        Assert.assertTrue(tmpHasBoth >= 0);
    }

    /**
     * This test method does some glycosylation statistics focused on linear sugars of COCONUT. It calculates:
     * How many molecules have linear sugars, only linear sugars (no circular ones), terminal linear sugars,
     * only terminal linear sugars (no non-terminal linear ones), non-terminal linear sugars, only non-terminal
     * linear sugars (no terminal linear ones), terminal and non-terminal linear sugars, and linear sugars in rings.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_linear_sugars_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsLinearSugarsTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_linear_sugars_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasNonTerminalLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasOnlyTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasOnlyNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasTerminalAndNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalLinearSugarsCNPs = new ArrayList<>(10);
        int tmpHasOnlyLinearSugarsCounter = 0;
        List<String> tmpHasOnlyLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasLinearSugarsInRingCounter = 0;
        List<String> tmpHasLinearSugarsInRingCNPs = new ArrayList<>(500);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpSugarRemovalUtil.hasLinearSugars(tmpMolecule);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpSugarRemovalUtil.hasCircularSugars(tmpMolecule);
                if (tmpHasAnyLinearSugar) {
                    tmpHasLinearSugarsCounter++;
                    tmpHasLinearSugarsCNPs.add(tmpID);
                    if (!tmpHasAnyCircularSugar) {
                        tmpHasOnlyLinearSugarsCounter++;
                        tmpHasOnlyLinearSugarsCNPs.add(tmpID);
                    }
                    //terminal and non-terminal
                    List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());
                    int tmpNumberOfLinearSugarMoieties;
                    int tmpNumberOfTerminalLinearSugarMoieties;
                    int tmpNumberOfNonTerminalLinearSugarMoieties;
                    tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();
                    //note:linear moieties that become terminal after removal of a circular moiety are not counted here!
                    List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                    //-1 for the deglycosylated core at the beginning of the list
                    tmpNumberOfTerminalLinearSugarMoieties = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                    tmpNumberOfNonTerminalLinearSugarMoieties = tmpNumberOfLinearSugarMoieties - tmpNumberOfTerminalLinearSugarMoieties;
                    Assert.assertTrue(tmpNumberOfNonTerminalLinearSugarMoieties >= 0);
                    if (tmpNumberOfTerminalLinearSugarMoieties > 0) {
                        tmpHasTerminalLinearSugarsCounter++;
                        tmpHasTerminalLinearSugarsCNPs.add(tmpID);
                        if (tmpNumberOfNonTerminalLinearSugarMoieties == 0) {
                            tmpHasOnlyTerminalLinearSugarsCounter++;
                            tmpHasOnlyTerminalLinearSugarsCNPs.add(tmpID);
                        }
                    }
                    if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                        tmpHasNonTerminalLinearSugarsCounter++;
                        tmpHasNonTerminalLinearSugarsCNPs.add(tmpID);
                        if (tmpNumberOfTerminalLinearSugarMoieties == 0) {
                            tmpHasOnlyNonTerminalLinearSugarsCounter++;
                            tmpHasOnlyNonTerminalLinearSugarsCNPs.add(tmpID);
                        }
                    }
                    if (tmpNumberOfTerminalLinearSugarMoieties > 0 && tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                        tmpHasTerminalAndNonTerminalLinearSugarsCounter++;
                        tmpHasTerminalAndNonTerminalLinearSugarsCNPs.add(tmpID);
                    }
                }
                //leaving default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                IAtomContainer tmpNewClone = tmpMolecule.clone();
                List<IAtomContainer> tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                int tmpListSizeWithCandidatesInCycles = tmpLinearCandidates.size();
                if (tmpListSizeWithCandidatesInCycles > 0) {
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    int tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidates.size();
                    int tmpNumberOfLinearSugarsInCycles = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    if (tmpNumberOfLinearSugarsInCycles > 0) {
                        tmpHasLinearSugarsInRingCounter++;
                        tmpHasLinearSugarsInRingCNPs.add(tmpID);
                    }
                }
                //back to default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        double tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Only linear sugar containing molecules counter: " + tmpHasOnlyLinearSugarsCounter);
        tmpOutputWriter.println("Only linear sugar containing molecules counter: " + tmpHasOnlyLinearSugarsCounter);
        System.out.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        System.out.println("Only terminal linear sugar containing molecules counter: " + tmpHasOnlyTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Only terminal linear sugar containing molecules counter: " + tmpHasOnlyTerminalLinearSugarsCounter);
        System.out.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        System.out.println("Only non-terminal linear sugar containing molecules counter: " + tmpHasOnlyNonTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Only non-terminal linear sugar containing molecules counter: " + tmpHasOnlyNonTerminalLinearSugarsCounter);
        System.out.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        System.out.println("Linear sugar moieties in rings containing molecules counter: " + tmpHasLinearSugarsInRingCounter);
        tmpOutputWriter.println("Linear sugar moieties in rings containing molecules counter: " + tmpHasLinearSugarsInRingCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Linear sugar containing molecules: " + tmpHasLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only linear sugar containing molecules: " + tmpHasOnlyLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal linear sugars containing molecules: " + tmpHasTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only terminal linear sugar containing molecules: " + tmpHasOnlyTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Non-terminal linear sugar containing molecules: " + tmpHasNonTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Only non-terminal linear sugar containing molecules: " + tmpHasOnlyNonTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal and non-terminal linear sugar containing molecules: " + tmpHasTerminalAndNonTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Linear sugar moieties in rings containing molecules: " + tmpHasLinearSugarsInRingCNPs);
        tmpOutputWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpHasTerminalAndNonTerminalLinearSugarsCounter
                + tmpHasOnlyTerminalLinearSugarsCounter
                + tmpHasOnlyNonTerminalLinearSugarsCounter);
    }

    /**
     * This test method does some glycosylation statistics focused on circular sugar moieties detected in COCONUT.
     * It calculates: How many circular sugar moieties, terminal circular sugar moieties, non-terminal circular sugar
     * moieties, circular sugar moieties with glycosidic bonds, and terminal/non-terminal circular sugar moieties with
     * glycosidic bonds can be detected in COCONUT. Also, the numbers of detected furanoses, pyranoses, and heptoses
     * are reported and compiled in a CSV file.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_circular_sugar_moieties_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsCircularSugarMoietiesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_circular_sugar_moieties_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVWriter = this.initializeOutputFile(tmpOutputFolderPath, "CircSugarsSizeFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpCircularSugarMoietiesCounter = 0;
        int tmpTerminalCircularSugarMoietiesCounter = 0;
        int tmpNonTerminalCircularSugarMoietiesCounter = 0;
        int tmpCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest have no glycosidic bond
        int tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest are non-terminal
        HashMap<Integer, Integer> tmpFrequenciesOfSizesOfCircularSugarMoietiesMap = new HashMap<>(10, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpSugarRemovalUtil.hasCircularSugars(tmpMolecule);
                if (tmpHasAnyCircularSugar) {
                    //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                    List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                    int tmpNumberOfCircularSugarMoieties;
                    int tmpNumberOfTerminalCircularSugarMoieties;
                    int tmpNumberOfNonTerminalCircularSugarMoieties;
                    int tmpNumberOfGlycosidicBonds;
                    int tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                    int tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond;
                    tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();
                    //the first is the overall counter, the second one is specific for this molecule
                    tmpCircularSugarMoietiesCounter += tmpNumberOfCircularSugarMoieties;
                    for (IAtomContainer tmpCircularSugarCandidate : tmpCircularSugarCandidatesList) {
                        int tmpCandidateSize = AtomContainerManipulator.getHeavyAtoms(tmpCircularSugarCandidate).size();
                        if (!tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.containsKey(tmpCandidateSize)) {
                            tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.put(tmpCandidateSize, 1);
                        } else {
                            Integer tmpCurrentCount = tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(tmpCandidateSize);
                            tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.put(tmpCandidateSize, tmpCurrentCount + 1);
                        }
                    }
                    //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                    List<IAtomContainer> tmpRemovedTerminalCircularSugarMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                    //-1 for the deglycosylated core at the beginning of the list
                    tmpNumberOfTerminalCircularSugarMoieties = tmpRemovedTerminalCircularSugarMoieties.size() - 1 ;
                    tmpNumberOfNonTerminalCircularSugarMoieties = tmpNumberOfCircularSugarMoieties - tmpNumberOfTerminalCircularSugarMoieties;
                    Assert.assertTrue(tmpNumberOfNonTerminalCircularSugarMoieties >= 0);
                    //leaving default! Now, only circular sugars having glycosidic bonds are in the candidates and removed moieties
                    tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
                    boolean tmpMoleculeQualifiesForExemption = tmpSugarRemovalUtil.isQualifiedForGlycosidicBondExemption(tmpMolecule.clone());
                    if (tmpMoleculeQualifiesForExemption) {
                        tmpNumberOfGlycosidicBonds = 0;
                        tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = 0;
                    } else {
                        List<IAtomContainer> tmpCircularSugarCandidatesWithGlycosidicBondsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                        tmpNumberOfGlycosidicBonds = tmpCircularSugarCandidatesWithGlycosidicBondsList.size();
                        List<IAtomContainer> tmpRemovedTerminalCircularMoietiesWithGlycosidicBond = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond = tmpRemovedTerminalCircularMoietiesWithGlycosidicBond.size() - 1;
                        tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond = tmpNumberOfGlycosidicBonds - tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                        Assert.assertTrue(tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond >= 0);
                    }
                    //back to default!
                    tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
                    if (tmpNumberOfTerminalCircularSugarMoieties > 0) {
                        tmpTerminalCircularSugarMoietiesCounter += tmpNumberOfTerminalCircularSugarMoieties;
                    }
                    if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                        tmpNonTerminalCircularSugarMoietiesCounter += tmpNumberOfNonTerminalCircularSugarMoieties;
                    }
                    if (tmpNumberOfGlycosidicBonds > 0) {
                        //the first is the overall counter, the second one is specific for this molecule
                        tmpCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfGlycosidicBonds;
                        if (tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond > 0) {
                            //the first is the overall counter, the second one is specific for this molecule
                            tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                        }
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Detected circular sugar moieties counter: " + tmpCircularSugarMoietiesCounter);
        tmpOutputWriter.println("Detected circular sugar moieties counter: " + tmpCircularSugarMoietiesCounter);
        System.out.println("Detected terminal circular sugar moieties counter: " + tmpTerminalCircularSugarMoietiesCounter);
        tmpOutputWriter.println("Detected terminal circular sugar moieties counter: " + tmpTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected non-terminal circular sugar moieties counter: " + tmpNonTerminalCircularSugarMoietiesCounter);
        tmpOutputWriter.println("Detected non-terminal circular sugar moieties counter: " + tmpNonTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected circular sugar moieties that have a glycosidic bond counter: " + tmpCircularSugarMoietiesWithGlycosidicBondCounter);
        tmpOutputWriter.println("Detected circular sugar moieties that have a glycosidic bond counter: " + tmpCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesCounter - tmpCircularSugarMoietiesWithGlycosidicBondCounter) + " circular sugar moieties do not have a glycosidic bond.");
        tmpOutputWriter.println((tmpCircularSugarMoietiesCounter - tmpCircularSugarMoietiesWithGlycosidicBondCounter) + " circular sugar moieties do not have a glycosidic bond.");
        System.out.println("Detected circular sugar moieties that have a glycosidic bond and are terminal counter: " + tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter);
        tmpOutputWriter.println("Detected circular sugar moieties that have a glycosidic bond and are terminal counter: " + tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesWithGlycosidicBondCounter - tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter)
                + " circular sugar moieties that have a glycosidic bond are non-terminal.");
        tmpOutputWriter.println((tmpCircularSugarMoietiesWithGlycosidicBondCounter - tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter)
                + " circular sugar moieties that have a glycosidic bond are non-terminal.");
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Size (= heavy atom count) frequency distribution of circular sugars: ");
        tmpOutputWriter.println("Size (= heavy atom count) frequency distribution of circular sugars: ");
        tmpCSVWriter.println("HeavyAtomCount" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfCircularSugars = 0;
        int tmpMaxHeavyAtomCount = Collections.max(tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.keySet());
        //starting at zero to see whether problems occurred
        for (int i = 0; i <= tmpMaxHeavyAtomCount; i++) {
            Integer tmpFrequency = tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfCircularSugars += tmpFrequency;
        }
        tmpOutputWriter.flush();
        tmpCSVWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVWriter.close();
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTerminalCircularSugarMoietiesCounter
                + tmpNonTerminalCircularSugarMoietiesCounter);
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTotalOfCircularSugars);
    }

    /**
     * This test method does some glycosylation statistics focused on linear sugar moieties detected in COCONUT.
     * It calculates: How many linear sugar moieties, terminal linear sugar moieties, and non-terminal linear sugar
     * moieties can be detected in COCONUT. Also, the numbers of detected tetroses, pentoses, hexoses and heptoses
     * are reported and compiled in a CSV file.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_linear_sugar_moieties_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsLinearSugarMoietiesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_linear_sugar_moieties_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVheavyAtomCountWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugarsHeavyAtomCountFrequencies.csv");
        PrintWriter tmpCSVcarbonAtomCountWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugarsCarbonAtomCountFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpLinearSugarMoietiesCounter = 0;
        int tmpTerminalLinearSugarMoietiesCounter = 0;
        int tmpNonTerminalLinearSugarMoietiesCounter = 0;
        HashMap<Integer, Integer> tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap = new HashMap<>(10, 0.9f);
        int tmpLinSugLostInRemovalOfCircSugCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpSugarRemovalUtil.hasLinearSugars(tmpMolecule);
                if (tmpHasAnyLinearSugar) {
                    //terminal and non-terminal
                    List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());
                    int tmpNumberOfLinearSugarMoieties;
                    int tmpNumberOfTerminalLinearSugarMoieties;
                    int tmpNumberOfNonTerminalLinearSugarMoieties;
                    tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();
                    //the first is the overall counter, the second one is specific for this molecule
                    tmpLinearSugarMoietiesCounter += tmpNumberOfLinearSugarMoieties;
                    for (IAtomContainer tmpLinearSugarCandidate : tmpLinearSugarCandidatesList) {
                        int tmpCandidateSize = AtomContainerManipulator.getHeavyAtoms(tmpLinearSugarCandidate).size();
                        if (!tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.containsKey(tmpCandidateSize)) {
                            tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.put(tmpCandidateSize, 1);
                        } else {
                            Integer tmpCurrentCount = tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(tmpCandidateSize);
                            tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.put(tmpCandidateSize, tmpCurrentCount + 1);
                        }
                        int tmpCarbonCount = 0;
                        for (IAtom tmpAtom : tmpLinearSugarCandidate.atoms()) {
                            String tmpSymbol = tmpAtom.getSymbol();
                            if (tmpSymbol.equals("C")) {
                                tmpCarbonCount++;
                            }
                        }
                        if (!tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.containsKey(tmpCarbonCount)) {
                            tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, 1);
                        } else {
                            Integer tmpCurrentCount = tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(tmpCarbonCount);
                            tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, tmpCurrentCount + 1);
                        }
                    }
                    //note:linear moieties that become terminal after removal of a circular moiety are not counted here!
                    List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                    //-1 for the deglycosylated core at the beginning of the list
                    tmpNumberOfTerminalLinearSugarMoieties = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                    tmpNumberOfNonTerminalLinearSugarMoieties = tmpNumberOfLinearSugarMoieties - tmpNumberOfTerminalLinearSugarMoieties;
                    Assert.assertTrue(tmpNumberOfNonTerminalLinearSugarMoieties >= 0);
                    if (tmpNumberOfTerminalLinearSugarMoieties > 0) {
                        tmpTerminalLinearSugarMoietiesCounter += tmpNumberOfTerminalLinearSugarMoieties;
                    }
                    if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                        tmpNonTerminalLinearSugarMoietiesCounter += tmpNumberOfNonTerminalLinearSugarMoieties;
                    }
                    //Testing for linear sugar candidates that get lost through the removal of circular sugars;
                    // NOT including those in circles, see below
                    IAtomContainer tmpNewClone = tmpMolecule.clone();
                    //leaving default!
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
                    tmpSugarRemovalUtil.removeCircularSugars(tmpNewClone, false);
                    List<IAtomContainer> tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    int tmpLinSugLostInRemovalOfCircSug = tmpNumberOfLinearSugarMoieties - tmpLinearCandidates.size();
                    if (tmpLinSugLostInRemovalOfCircSug < 0) {
                        throw new Exception();
                    }
                    tmpLinSugLostInRemovalOfCircSugCounter += tmpLinSugLostInRemovalOfCircSug;
                    //restoring default!
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Detected linear sugar moieties counter: " + tmpLinearSugarMoietiesCounter);
        tmpOutputWriter.println("Detected linear sugar moieties counter: " + tmpLinearSugarMoietiesCounter);
        System.out.println("Detected terminal linear sugar moieties counter: " + tmpTerminalLinearSugarMoietiesCounter);
        tmpOutputWriter.println("Detected terminal linear sugar moieties counter: " + tmpTerminalLinearSugarMoietiesCounter);
        System.out.println("Detected non-terminal linear sugar moieties counter: " + tmpNonTerminalLinearSugarMoietiesCounter);
        tmpOutputWriter.println("Detected non-terminal linear sugar moieties counter: " + tmpNonTerminalLinearSugarMoietiesCounter);
        System.out.println("Number of detected linear sugars that got lost through " +
                "the removal of circular sugars counter: " + tmpLinSugLostInRemovalOfCircSugCounter);
        tmpOutputWriter.println("Number of detected linear sugars that got lost through " +
                "the removal of circular sugars counter: " + tmpLinSugLostInRemovalOfCircSugCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Size (= heavy atom count) frequency distribution of linear sugars: ");
        tmpOutputWriter.println("Size (= heavy atom count) frequency distribution of linear sugars: ");
        tmpCSVheavyAtomCountWriter.println("HeavyAtomCount" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfLinearSugars1 = 0;
        int tmpMaxHeavyAtomCount = Collections.max(tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.keySet());
        //starting at zero to see whether problems occurred
        for (int i = 0; i <= tmpMaxHeavyAtomCount; i++) {
            Integer tmpFrequency = tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + ":" + tmpFrequency);
            tmpOutputWriter.println(i + ":" + tmpFrequency);
            tmpCSVheavyAtomCountWriter.println(i + ":" + tmpFrequency);
            tmpTotalOfLinearSugars1 += tmpFrequency;
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Size (= carbon atom count) frequency distribution of linear sugars (note set min and max sizes): ");
        tmpOutputWriter.println("Size (= carbon atom count) frequency distribution of linear sugars (note set min and max sizes): ");
        tmpCSVcarbonAtomCountWriter.println("CarbonAtomCount" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfLinearSugars2 = 0;
        int tmpMaxCarbonAtomCount = Collections.max(tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.keySet());
        //starting at zero to see whether problems occurred
        for (int i = 0; i <= tmpMaxCarbonAtomCount; i++) {
            Integer tmpFrequency = tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVcarbonAtomCountWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfLinearSugars2 += tmpFrequency;
        }
        tmpOutputWriter.flush();
        tmpCSVheavyAtomCountWriter.flush();
        tmpCSVcarbonAtomCountWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVheavyAtomCountWriter.close();
        tmpCSVcarbonAtomCountWriter.close();
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTerminalLinearSugarMoietiesCounter
                + tmpNonTerminalLinearSugarMoietiesCounter);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars1);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars2);
    }

    /**
     * This test method does some glycosylation statistics focused on linear sugar moieties that are ring substructures
     * detected in COCONUT.
     * It calculates how many linear sugar moieties that are substructures of cycles can be detected in COCONUT.
     * Also, their size distribution is reported and compiled in a CSV file. Additionally, images of every atom with such
     * cyclic linear sugar moieties are created in the output directory where the moieties are highlighted.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_linear_sugar_moieties_in_rings_test.
     * To run this test, remove the JUnit @Ignore tag that is set because
     * this test creates many files.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void coconutStatsLinearSugarMoietiesInRingsTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_linear_sugar_moieties_in_rings_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugarsCarbonAtomCountFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpLinearSugarMoietiesInRingsCounter = 0;
        int tmpLinSugInRingsLostInRemovalOfCircSugCounter = 0;
        List<String> tmpLinSugInRingsLostInRemovalOfCircSugCNPs = new ArrayList<>(60);
        HashMap<Integer, Integer> tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap = new HashMap<>(10, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //leaving default SRU settings
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                //note: per default, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpSugarRemovalUtil.hasLinearSugars(tmpMolecule);
                if (tmpHasAnyLinearSugar) {
                    List<IAtomContainer> tmpLinearCandidatesIncludingCycles = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule);
                    int tmpListSizeWithCandidatesInCycles = tmpLinearCandidatesIncludingCycles.size();
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    List<IAtomContainer> tmpLinearCandidatesExcludingCycles = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule);
                    int tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidatesExcludingCycles.size();
                    int tmpNumberOfLinearSugarsInCycles = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    Assert.assertTrue(tmpNumberOfLinearSugarsInCycles >= 0);
                    if (tmpNumberOfLinearSugarsInCycles > 0) {
                        tmpLinearSugarMoietiesInRingsCounter += tmpNumberOfLinearSugarsInCycles;
                        //this will be the list of candidates that only get detected if cyclic atoms are included
                        List<IAtomContainer> tmpLinearCandidatesActuallyInRings = new ArrayList<>(tmpNumberOfLinearSugarsInCycles * 2);
                        int[][] tmpAdjList = GraphUtil.toAdjList(tmpMolecule);
                        RingSearch tmpRingSearch = new RingSearch(tmpMolecule, tmpAdjList);
                        for (IAtomContainer tmpCandidate : tmpLinearCandidatesIncludingCycles) {
                            for (IAtom tmpAtom : tmpCandidate.atoms()) {
                                if (tmpRingSearch.cyclic(tmpAtom)) {
                                    int tmpCarbonCount = 0;
                                    for (IAtom tmpAtom2 : tmpCandidate.atoms()) {
                                        String tmpSymbol = tmpAtom2.getSymbol();
                                        if (tmpSymbol.equals("C")) {
                                            tmpCarbonCount++;
                                        }
                                    }
                                    if (!tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.containsKey(tmpCarbonCount)) {
                                        tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, 1);
                                    } else {
                                        Integer tmpCurrentCount = tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(tmpCarbonCount);
                                        tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.put(tmpCarbonCount, tmpCurrentCount + 1);
                                    }
                                    tmpLinearCandidatesActuallyInRings.add(tmpCandidate);
                                    //move on with the next candidate, i.e. break the loop over the atoms in this candidate
                                    break;
                                }
                            }
                        }
                        tmpDepictionGenerator.withHighlight(tmpLinearCandidatesActuallyInRings, Color.BLUE)
                                .withSize(2000, 2000)
                                .withFillToFit()
                                .depict(tmpMolecule)
                                .writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    }
                    //leaving default further!
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
                    IAtomContainer tmpNewClone = tmpMolecule.clone();
                    tmpSugarRemovalUtil.removeCircularSugars(tmpNewClone, false);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    tmpLinearCandidatesIncludingCycles = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpListSizeWithCandidatesInCycles = tmpLinearCandidatesIncludingCycles.size();
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidatesExcludingCycles = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidatesExcludingCycles.size();
                    int tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    int tmpLinSugInRingsLostInRemovalOfCircSug = tmpNumberOfLinearSugarsInCycles - tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars;
                    Assert.assertTrue(tmpLinSugInRingsLostInRemovalOfCircSug >= 0);
                    if (tmpLinSugInRingsLostInRemovalOfCircSug > 0) {
                        tmpLinSugInRingsLostInRemovalOfCircSugCounter += tmpLinSugInRingsLostInRemovalOfCircSug;
                        tmpLinSugInRingsLostInRemovalOfCircSugCNPs.add(tmpID);
                    }
                    //back to this default
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
                }
                //back to default settings
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        //Note: these are moieties that were not detected while excluding circular atoms, so i.e. they are 'mainly' in rings
        // other linear moieties might include circular atoms but be big enough to be detected without them
        System.out.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter);
        tmpOutputWriter.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter);
        System.out.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinSugInRingsLostInRemovalOfCircSugCounter);
        tmpOutputWriter.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinSugInRingsLostInRemovalOfCircSugCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules that lost a linear sugar in a ring after removal of circular sugars: "
                + tmpLinSugInRingsLostInRemovalOfCircSugCNPs);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Size (= carbon atom count) frequency distribution of linear sugars (note set min and max sizes): ");
        tmpOutputWriter.println("Size (= carbon atom count) frequency distribution of linear sugars (note set min and max sizes): ");
        tmpCSVWriter.println("CarbonAtomCount" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpMaxCarbonAtomCount = Collections.max(tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.keySet());
        //starting at zero to see whether problems occurred
        for (int i = 0; i <= tmpMaxCarbonAtomCount; i++) {
            Integer tmpFrequency = tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
        }
        tmpOutputWriter.flush();
        tmpCSVWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVWriter.close();
    }

    /**
     * This test method does some glycosylation statistics focused on molecules that are sugar monomers or polymers in
     * COCONUT. It calculates how many molecules are basically sugars (i.e. consist only of sugar units) and what
     * portions of these are circular or linear sugar polymers or monomers. All statistics
     * are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_sugar_molecules_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsSugarMoleculesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_sugar_molecules_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        //number of molecules that are basically sugars, circular or linear, polymer or single unit
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarCNPs = new ArrayList<>(2000);
        int tmpBasicallyACircularSugarCounter = 0;
        List<String> tmpBasicallyACircularSugarCNPs = new ArrayList<>(2000);
        int tmpBasicallyALinearSugarCounter = 0;
        List<String> tmpBasicallyALinearSugarCNPs = new ArrayList<>(2000);
        int tmpBasicallyASingleSugarUnitCounter = 0; //circular or linear
        List<String> tmpBasicallyASingleSugarUnitCNPs = new ArrayList<>(1000);
        int tmpBasicallyASingleCircularSugarCounter = 0;
        List<String> tmpBasicallyASingleCircularSugarCNPs = new ArrayList<>(1000);
        int tmpBasicallyASingleLinearSugarCounter = 0;
        List<String> tmpBasicallyASingleLinearSugarCNPs = new ArrayList<>(1000);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //removes only terminal moieties but that is correct here
                List<IAtomContainer> tmpDeglycosylatedCloneAndRemovedSugarMoietiesList =
                        tmpSugarRemovalUtil.removeAndReturnCircularAndLinearSugars(tmpMolecule, true);
                IAtomContainer tmpDeglycosylatedClone = tmpDeglycosylatedCloneAndRemovedSugarMoietiesList.get(0);
                if (tmpDeglycosylatedClone.isEmpty()) {
                    tmpBasicallyASugarCounter++;
                    tmpBasicallyASugarCNPs.add(tmpID);
                    //note: it is important to count the actually removed moieties here, not the detected ones!
                    // Because there are multiple rounds of detection in the removal if only terminal moieties are removed
                    int tmpNumberOfMoieties = tmpDeglycosylatedCloneAndRemovedSugarMoietiesList.size() - 1;
                    if (tmpNumberOfMoieties == 1) {
                        tmpBasicallyASingleSugarUnitCounter++;
                        tmpBasicallyASingleSugarUnitCNPs.add(tmpID);
                    }
                }
                IAtomContainer tmpCircularDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularSugars(tmpMolecule, true);
                if (tmpCircularDeglycosylatedClone.isEmpty()) {
                    tmpBasicallyACircularSugarCounter++;
                    tmpBasicallyACircularSugarCNPs.add(tmpID);
                    //note: here, it is ok to only count the detected moieties because there is only one round of detection in the removal
                    int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfCircularSugars(tmpMolecule.clone());
                    if (tmpNumberOfMoieties == 1) {
                        tmpBasicallyASingleCircularSugarCounter++;
                        tmpBasicallyASingleCircularSugarCNPs.add(tmpID);
                    }
                }
                IAtomContainer tmpLinearDeglycosylatedClone = tmpSugarRemovalUtil.removeLinearSugars(tmpMolecule, true);
                if (tmpLinearDeglycosylatedClone.isEmpty()) {
                    tmpBasicallyALinearSugarCounter++;
                    tmpBasicallyALinearSugarCNPs.add(tmpID);
                    //note: here, it is ok to only count the detected moieties because there is only one round of detection in the removal
                    int tmpNumberOfMoieties = tmpSugarRemovalUtil.getNumberOfLinearSugars(tmpMolecule);
                    if (tmpNumberOfMoieties == 1) {
                        tmpBasicallyASingleLinearSugarCounter++;
                        tmpBasicallyASingleLinearSugarCNPs.add(tmpID);
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        System.out.println("Basically a single sugar unit counter: " + tmpBasicallyASingleSugarUnitCounter);
        tmpOutputWriter.println("Basically a single sugar unit counter: " + tmpBasicallyASingleSugarUnitCounter);
        System.out.println("Basically a circular sugar counter: " + tmpBasicallyACircularSugarCounter);
        tmpOutputWriter.println("Basically a circular sugar counter: " + tmpBasicallyACircularSugarCounter);
        System.out.println("Basically a single circular sugar counter: " + tmpBasicallyASingleCircularSugarCounter);
        tmpOutputWriter.println("Basically a single circular sugar counter: " + tmpBasicallyASingleCircularSugarCounter);
        System.out.println("Basically a linear sugar counter: " + tmpBasicallyALinearSugarCounter);
        tmpOutputWriter.println("Basically a linear sugar counter: " + tmpBasicallyALinearSugarCounter);
        System.out.println("Basically a single linear sugar counter: " + tmpBasicallyASingleLinearSugarCounter);
        tmpOutputWriter.println("Basically a single linear sugar counter: " + tmpBasicallyASingleLinearSugarCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a sugar: " + tmpBasicallyASugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a single sugar unit: " + tmpBasicallyASingleSugarUnitCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a circular sugar: " + tmpBasicallyACircularSugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a single circular sugar: " + tmpBasicallyASingleCircularSugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a linear sugar: " + tmpBasicallyALinearSugarCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a single linear sugar: " + tmpBasicallyASingleLinearSugarCNPs);
        tmpOutputWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpBasicallyASingleSugarUnitCounter,
                tmpBasicallyASingleCircularSugarCounter + tmpBasicallyASingleLinearSugarCounter);
    }


    /**
     * This test method does some glycosylation statistics focused on the general distribution of sugar moieties detected
     * in COCONUT. It calculates how many molecules have how many (circular/linear) sugar moieties, respectively.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_moiety_frequencies_test. Also, CSV files for the sugar moiety
     * distributions are created.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsSugarMoietiesFrequenciesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_moiety_frequencies_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVmoietyNrFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "SugMoietyNrFrequencies.csv");
        PrintWriter tmpCSVcircMoietyNrFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "CircSugMoietyNrFrequencies.csv");
        PrintWriter tmpCSVcircMoietyGlyBondNrFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "CircSugMoietyGlyBondNrFrequencies.csv");
        PrintWriter tmpCSVlinMoietyNrFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugMoietyNrFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        int tmpHasGlycosidicBondCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManySugarMoietiesMap = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap = new HashMap<>(10, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //using default settings where nothing else is specified
                boolean tmpHasAnyTypeOfSugar = tmpSugarRemovalUtil.hasCircularOrLinearSugars(tmpMolecule);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpHasAnyTypeOfSugarsCounter++;
                    int tmpNumberOfCircularAndLinearSugarMoieties = tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpMolecule);
                    if (!tmpHowManyMoleculesHaveHowManySugarMoietiesMap.containsKey(tmpNumberOfCircularAndLinearSugarMoieties)) {
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, 1);
                    } else {
                        Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(tmpNumberOfCircularAndLinearSugarMoieties);
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, tmpCurrentListValue + 1);
                    }
                    if (tmpHasAnyCircularSugar) {
                        tmpHasCircularSugarsCounter++;
                        //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                        List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfCircularSugarMoieties;
                        tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();
                        if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.containsKey(tmpNumberOfCircularSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(tmpNumberOfCircularSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, tmpCurrentValue + 1);
                        }
                        int tmpNumberOfGlycosidicBonds;
                        //leaving default! Now, only circular sugars having glycosidic bonds are in the candidates and removed moieties
                        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(true);
                        boolean tmpMoleculeQualifiesForExemption = tmpSugarRemovalUtil.isQualifiedForGlycosidicBondExemption(tmpMolecule.clone());
                        if (tmpMoleculeQualifiesForExemption) {
                            tmpNumberOfGlycosidicBonds = 0;
                        } else {
                            List<IAtomContainer> tmpCircularSugarCandidatesWithGlycosidicBondsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                            tmpNumberOfGlycosidicBonds = tmpCircularSugarCandidatesWithGlycosidicBondsList.size();
                        }
                        //back to default!
                        tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithOGlycosidicBondSetting(false);
                        if (tmpNumberOfGlycosidicBonds > 0) {
                            tmpHasGlycosidicBondCounter++;
                            if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.containsKey(tmpNumberOfGlycosidicBonds)) {
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, 1);
                            } else {
                                Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(tmpNumberOfGlycosidicBonds);
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, tmpCurrentListValue + 1);
                            }
                        }
                    }
                    if (tmpHasAnyLinearSugar) {
                        tmpHasLinearSugarsCounter++;
                        //terminal and non-terminal
                        List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfLinearSugarMoieties;
                        tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();
                        if (!tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.containsKey(tmpNumberOfLinearSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(tmpNumberOfLinearSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, tmpCurrentValue + 1);
                        }
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("How many molecules have how many sugars: ");
        tmpOutputWriter.println("How many molecules have how many sugars: ");
        tmpCSVmoietyNrFreqWriter.println("NrOfMoieties" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfSugarContainingMolecules = 0;
        int tmpMaxNrOfMoieties = Collections.max(tmpHowManyMoleculesHaveHowManySugarMoietiesMap.keySet());
        for (int i = 1; i <= tmpMaxNrOfMoieties; i++) {
            Integer tmpFrequency = tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVmoietyNrFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfSugarContainingMolecules += tmpFrequency;
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("How many molecules have how many circular sugars (out of these that have any): ");
        tmpOutputWriter.println("How many molecules have how many circular sugars (out of these that have any): ");
        tmpCSVcircMoietyNrFreqWriter.println("NrOfCircMoieties" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfCircularSugarContainingMolecules = 0;
        int tmpMaxNrOfCircMoieties = Collections.max(tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.keySet());
        for (int i = 1; i <= tmpMaxNrOfCircMoieties; i++) {
            Integer tmpFrequency = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVcircMoietyNrFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfCircularSugarContainingMolecules += tmpFrequency;
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("How many molecules have how many circular sugars attached via a glycosidic bond (out of these that have any): ");
        tmpOutputWriter.println("How many molecules have how many circular sugars attached via a glycosidic bond (out of these that have any): ");
        tmpCSVcircMoietyGlyBondNrFreqWriter.println("NrOfCircMoietiesWithGlyBond" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules = 0;
        int tmpMaxNrOfCircMoietiesGlyBond = Collections.max(tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.keySet());
        for (int i = 1; i <= tmpMaxNrOfCircMoietiesGlyBond; i++) {
            Integer tmpFrequency = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVcircMoietyGlyBondNrFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules += tmpFrequency;
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("How many molecules have how many linear sugars (out of these that have any): ");
        tmpOutputWriter.println("How many molecules have how many linear sugars (out of these that have any): ");
        tmpCSVlinMoietyNrFreqWriter.println("NrOfLinMoieties" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpTotalOfLinearSugarContainingMolecules = 0;
        int tmpMaxNrOfLinMoieties = Collections.max(tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.keySet());
        for (int i = 1; i <= tmpMaxNrOfLinMoieties; i++) {
            Integer tmpFrequency = tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVlinMoietyNrFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpTotalOfLinearSugarContainingMolecules += tmpFrequency;
        }
        tmpOutputWriter.flush();
        tmpCSVmoietyNrFreqWriter.flush();
        tmpCSVcircMoietyNrFreqWriter.flush();
        tmpCSVcircMoietyGlyBondNrFreqWriter.flush();
        tmpCSVlinMoietyNrFreqWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVmoietyNrFreqWriter.close();
        tmpCSVcircMoietyNrFreqWriter.close();
        tmpCSVcircMoietyGlyBondNrFreqWriter.close();
        tmpCSVlinMoietyNrFreqWriter.close();
        Assert.assertEquals(tmpHasAnyTypeOfSugarsCounter, tmpTotalOfSugarContainingMolecules);
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpTotalOfCircularSugarContainingMolecules);
        Assert.assertEquals(tmpHasGlycosidicBondCounter, tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules);
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpTotalOfLinearSugarContainingMolecules);
    }

    /**
     * This test method does some glycosylation statistics focused on the numbers and ratios of exocyclic oxygen atoms
     * of sugar moieties detected in COCONUT. It calculates the overall distribution of the exocyclic oxygen ratio (a
     * setting of the SRU) values of all types of circular sugar moieties and the distributions of exocyclic oxygen atom
     * counts for every type of circular sugars (furanoses, pyranoses, heptoses), respectively.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_exocyclic_oxygens_test. Also, CSV files for all
     * statistics are created.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsExocyclicOxygensTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_exocyclic_oxygens_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVExoCycOxFuranosesFreqWriter = this.initializeOutputFile(tmpOutputFolderPath,
                "ExoCycOxFuranosesFrequencies.csv");
        PrintWriter tmpCSVExoCycOxPyranosesFreqWriter = this.initializeOutputFile(tmpOutputFolderPath,
                "ExoCycOxPyranosesFrequencies.csv");
        PrintWriter tmpCSVExoCycOxHeptosesFreqWriter = this.initializeOutputFile(tmpOutputFolderPath,
                "ExoCycOxHeptosesFrequencies.csv");
        PrintWriter tmpCSVExoCycOxRatioFreqWriter = this.initializeOutputFile(tmpOutputFolderPath,
                "ExoCycOxRatioFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        NumberFormat tmpRatioOutputFormat = NumberFormat.getInstance(Locale.US);
        tmpRatioOutputFormat.setMaximumFractionDigits(1);
        tmpRatioOutputFormat.setRoundingMode(RoundingMode.DOWN);
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap
                = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap
                = new HashMap<>(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap
                = new HashMap<>(10, 0.9f);
        HashMap<String, Integer> tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap
                = new HashMap<>(10, 0.9f);
        int tmpUnexpectedRingSizeCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //leaving default settings!
                tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(false);
                IAtomContainer tmpMoleculeClone = tmpMolecule.clone();
                List<IAtomContainer> tmpCircularSugarsList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMoleculeClone);
                if (tmpCircularSugarsList.size() > 0) {
                    for (IAtomContainer tmpCandidate : tmpCircularSugarsList) {
                        int tmpRingSize = tmpCandidate.getAtomCount();
                        int tmpExocyclicOxygenAtomsCount = this.getExocyclicOxygenAtomCount(tmpCandidate, tmpMoleculeClone);
                        double tmpAttachedOxygensToAtomsInRingRatio =
                                ((double) tmpExocyclicOxygenAtomsCount / (double) tmpRingSize);
                        //note: the ratios are not rounded, the remaining decimals are neglected (rounding mode 'down'),
                        // which is correct here because the respective setting is a threshold
                        String tmpRoundedRatio = tmpRatioOutputFormat.format(tmpAttachedOxygensToAtomsInRingRatio);
                        if (!tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.containsKey(tmpRoundedRatio)) {
                            tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.put(tmpRoundedRatio, 1);
                        } else {
                            Integer tmpCurrentCount = tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(tmpRoundedRatio);
                            tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.put(tmpRoundedRatio, tmpCurrentCount + 1);
                        }
                        switch (tmpRingSize) {
                            case 5:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            case 6:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            case 7:
                                if (!tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.containsKey(tmpExocyclicOxygenAtomsCount)) {
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, 1);
                                } else {
                                    Integer tmpCurrentCount = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(tmpExocyclicOxygenAtomsCount);
                                    tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.put(tmpExocyclicOxygenAtomsCount, tmpCurrentCount + 1);
                                }
                                break;
                            default:
                                tmpUnexpectedRingSizeCounter++;
                                break;
                        }
                    }
                }
                //back to default settings!
                tmpSugarRemovalUtil.setDetectCircularSugarsOnlyWithEnoughExocyclicOxygenAtomsSetting(true);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Frequency distribution of exocyclic oxygen atoms to atoms in ring ratios of circular sugars: ");
        tmpOutputWriter.println("Frequency distribution of exocyclic oxygen atoms to atoms in ring ratios of circular sugars: ");
        tmpCSVExoCycOxRatioFreqWriter.println("Ratio" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        //iterating using int and creating double in the loop is necessary because of pitfalls in double arithmetic
        double tmpMaxRatio = Double.parseDouble(Collections.max(tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.keySet()));
        int tmpMaxForLoop = Double.valueOf(tmpMaxRatio * 10).intValue();
        for (int i = 0; i <= tmpMaxForLoop; i++) {
            Integer tmpFrequency = tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(tmpRatioOutputFormat.format((double)i/10));
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            String tmpRatio = tmpRatioOutputFormat.format((double) i / 10);
            System.out.println(tmpRatio + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(tmpRatio + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVExoCycOxRatioFreqWriter.println(tmpRatio + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 5-membered circular sugars: ");
        tmpOutputWriter.println("Frequency distribution of exocyclic oxygen atom counts of 5-membered circular sugars: ");
        tmpCSVExoCycOxFuranosesFreqWriter.println("NrOfOxygens" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        int tmpMaxNrOfOxygens = Collections.max(tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.keySet());
        for (int i = 0; i <= tmpMaxNrOfOxygens; i++) {
            Integer tmpFrequency = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVExoCycOxFuranosesFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 6-membered circular sugars: ");
        tmpOutputWriter.println("Frequency distribution of exocyclic oxygen atom counts of 6-membered circular sugars: ");
        tmpCSVExoCycOxPyranosesFreqWriter.println("NrOfOxygens" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpMaxNrOfOxygens = Collections.max(tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.keySet());
        for (int i = 0; i <= tmpMaxNrOfOxygens; i++) {
            Integer tmpFrequency = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVExoCycOxPyranosesFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 7-membered circular sugars: ");
        tmpOutputWriter.println("Frequency distribution of exocyclic oxygen atom counts of 7-membered circular sugars: ");
        tmpCSVExoCycOxHeptosesFreqWriter.println("NrOfOxygens" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpMaxNrOfOxygens = Collections.max(tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.keySet());
        for (int i = 0; i <= tmpMaxNrOfOxygens; i++) {
            Integer tmpFrequency = tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i);
            if (Objects.isNull(tmpFrequency)) {
                tmpFrequency = 0;
            }
            System.out.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpOutputWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            tmpCSVExoCycOxHeptosesFreqWriter.println(i + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
        }
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Number of circular sugar moieties that had an unexpected ring size (should be zero!): " + tmpUnexpectedRingSizeCounter);
        tmpOutputWriter.println("Number of circular sugar moieties that had an unexpected ring size (should be zero!): " + tmpUnexpectedRingSizeCounter);
        tmpOutputWriter.flush();
        tmpCSVExoCycOxRatioFreqWriter.flush();
        tmpCSVExoCycOxFuranosesFreqWriter.flush();
        tmpCSVExoCycOxPyranosesFreqWriter.flush();
        tmpCSVExoCycOxHeptosesFreqWriter.flush();
        tmpOutputWriter.close();
        tmpCSVExoCycOxRatioFreqWriter.close();
        tmpCSVExoCycOxFuranosesFreqWriter.close();
        tmpCSVExoCycOxPyranosesFreqWriter.close();
        tmpCSVExoCycOxHeptosesFreqWriter.close();
    }

    /**
     * This test method does some glycosylation statistics focused on the non-terminal vs. terminal distinction of sugar
     * moieties detected in COCONUT.
     * It calculates: How many sugar moieties are terminal or non-terminal overall and also how many of these are
     * circular or linear, respectively.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_terminal_non-terminal_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsTerminalNonTerminalTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_terminal_non-terminal_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpCircularSugarMoietiesCounter = 0;
        int tmpLinearSugarMoietiesCounter = 0;
        int tmpTerminalCircularSugarMoietiesCounter = 0;
        int tmpNonTerminalCircularSugarMoietiesCounter = 0;
        int tmpTerminalLinearSugarMoietiesCounter = 0;
        int tmpNonTerminalLinearSugarMoietiesCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //using default settings where nothing else is specified
                boolean tmpHasAnyTypeOfSugar = tmpSugarRemovalUtil.hasCircularOrLinearSugars(tmpMolecule);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not counted!
                boolean tmpHasAnyLinearSugar = tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    if (tmpHasAnyCircularSugar) {
                        //terminal and non-terminal, having a glycosidic bond or not (see default settings)
                        List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfCircularSugarMoieties;
                        int tmpNumberOfTerminalCircularSugarMoieties;
                        int tmpNumberOfNonTerminalCircularSugarMoieties;
                        tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();
                        //the first is the overall counter, the second one is specific for this molecule
                        tmpCircularSugarMoietiesCounter += tmpNumberOfCircularSugarMoieties;
                        //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalCircularSugarMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        tmpNumberOfTerminalCircularSugarMoieties = tmpRemovedTerminalCircularSugarMoieties.size() - 1 ;
                        tmpNumberOfNonTerminalCircularSugarMoieties = tmpNumberOfCircularSugarMoieties - tmpNumberOfTerminalCircularSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalCircularSugarMoieties >= 0);
                        if (tmpNumberOfTerminalCircularSugarMoieties > 0) {
                            tmpTerminalCircularSugarMoietiesCounter += tmpNumberOfTerminalCircularSugarMoieties;
                        }
                        if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                            tmpNonTerminalCircularSugarMoietiesCounter += tmpNumberOfNonTerminalCircularSugarMoieties;
                        }
                    }
                    if (tmpHasAnyLinearSugar) {
                        //terminal and non-terminal
                        List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfLinearSugarMoieties;
                        int tmpNumberOfTerminalLinearSugarMoieties;
                        int tmpNumberOfNonTerminalLinearSugarMoieties;
                        tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();
                        //the first is the overall counter, the second one is specific for this molecule
                        tmpLinearSugarMoietiesCounter += tmpNumberOfLinearSugarMoieties;
                        //note: linear moieties that become terminal after removal of a circular moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        tmpNumberOfTerminalLinearSugarMoieties = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                        tmpNumberOfNonTerminalLinearSugarMoieties = tmpNumberOfLinearSugarMoieties - tmpNumberOfTerminalLinearSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalLinearSugarMoieties >= 0);
                        if (tmpNumberOfTerminalLinearSugarMoieties > 0) {
                            tmpTerminalLinearSugarMoietiesCounter += tmpNumberOfTerminalLinearSugarMoieties;
                        }
                        if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                            tmpNonTerminalLinearSugarMoietiesCounter += tmpNumberOfNonTerminalLinearSugarMoieties;
                        }
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Terminal Sugars counter: " + (tmpTerminalCircularSugarMoietiesCounter + tmpTerminalLinearSugarMoietiesCounter));
        tmpOutputWriter.println("Terminal Sugars counter: " + (tmpTerminalCircularSugarMoietiesCounter + tmpTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpTerminalCircularSugarMoietiesCounter + " of these are circular");
        tmpOutputWriter.println(tmpTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpTerminalLinearSugarMoietiesCounter + " of these are linear");
        tmpOutputWriter.println(tmpTerminalLinearSugarMoietiesCounter + " of these are linear");
        System.out.println("Non-terminal Sugars counter: " + (tmpNonTerminalCircularSugarMoietiesCounter + tmpNonTerminalLinearSugarMoietiesCounter));
        tmpOutputWriter.println("Non-terminal Sugars counter: " + (tmpNonTerminalCircularSugarMoietiesCounter + tmpNonTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpNonTerminalCircularSugarMoietiesCounter + " of these are circular");
        tmpOutputWriter.println(tmpNonTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpNonTerminalLinearSugarMoietiesCounter + " of these are linear");
        tmpOutputWriter.println(tmpNonTerminalLinearSugarMoietiesCounter + " of these are linear");
        tmpOutputWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTerminalCircularSugarMoietiesCounter
                + tmpNonTerminalCircularSugarMoietiesCounter);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTerminalLinearSugarMoietiesCounter
                + tmpNonTerminalLinearSugarMoietiesCounter);
    }

    /**
     * This test method performs substructure searching and match counting in COCONUT with the SRU linear sugar patterns used for
     * initial linear sugar moiety candidate detection.
     * As a result, it lists the respective pattern as a SMILES string and its occurrence in COCONUT molecules.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_linear_sugar_patterns_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsLinearSugarPatternsAppearanceTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_linear_sugar_patterns_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinearSugarPatternFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //to also test their appearance
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        //Note: Here, additional molecules could be added to the list to also test them
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugarPatternsList();
        List<List<Object>> tmpLinearSugarPatterns = new ArrayList<>(tmpLinearSugarsList.size());
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        for (String tmpLinearSugarString : tmpLinearSugarsList) {
            List<Object> tmpList = new ArrayList<>(4);
            tmpList.add(0, tmpLinearSugarString);
            tmpList.add(1, DfPattern.findSubstructure(tmpSmiPar.parseSmiles(tmpLinearSugarString)));
            tmpList.add(2, 0);
            tmpLinearSugarPatterns.add(tmpList);
        }
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                boolean tmpMolHasAMatch = false;
                for (List<Object> tmpEntry : tmpLinearSugarPatterns) {
                    DfPattern tmpPattern = (DfPattern) tmpEntry.get(1);
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpEntry.set(2, (int)tmpEntry.get(2) + 1);
                        tmpMolHasAMatch = true;
                    }
                }
                if (tmpMolHasAMatch) {
                    tmpHasLinearSugarsCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Molecules that have at least one match counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Molecules that have at least one match counter: " + tmpHasLinearSugarsCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Linear sugar pattern SMILES codes and detected frequencies: ");
        tmpOutputWriter.println("Linear sugar pattern SMILES codes and detected frequencies: ");
        tmpCSVWriter.println("SMILEScode" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        for (List<Object> tmpEntry : tmpLinearSugarPatterns) {
            System.out.println(tmpEntry.get(0) + " " + tmpEntry.get(2));
            tmpOutputWriter.println(tmpEntry.get(0) + " " + tmpEntry.get(2));
            tmpCSVWriter.println(tmpEntry.get(0) + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpEntry.get(2));
        }
        tmpOutputWriter.flush();
        tmpCSVWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVWriter.close();
    }

    /**
     * This test method compiles and counts the removed circular sugar moieties the SRU returns by removing sugars from
     * all COCONUT molecules. Image files of the individual moieties with their frequencies as names are created in the
     * directory ./GlycosylationStatisticsTest_Output/coconut_stats_removed_circular_moiety_frequencies_test. Note that
     * these numbers give the total appearance of the respective moiety, not how many molecules have this moiety.
     * Also note that the circular sugar moieties returned by the SRU only consist of the rings without any further
     * information/side-chains/substituents. Therefore, this test basically compiles the frequencies of furanoses,
     * pyranoses, and heptoses that get removed from COCONUT molecules by the SRU.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_removed_circular_moiety_frequencies_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsRemovedCircularMoietyFrequenciesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_removed_circular_moiety_frequencies_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVmoietyFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "CircSugMoietyFrequencies.csv");
        PrintWriter tmpCSVperMoleculeWriter = this.initializeOutputFile(tmpOutputFolderPath, "MoleculesWithCircSugars.csv");
        String tmpCSVperMoleculeFileHeader = "ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "originalMoleculeSMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "SugarMoietySMILES";
        tmpCSVperMoleculeWriter.println(tmpCSVperMoleculeFileHeader);
        tmpCSVperMoleculeWriter.flush();
        String tmpCSVmoietyFreqFileHeader = "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "frequency" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "firstOrigin";
        tmpCSVmoietyFreqWriter.println(tmpCSVmoietyFreqFileHeader);
        tmpCSVmoietyFreqWriter.flush();
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        String tmpOutput = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        int tmpDifferentMoietiesCounter = 0;
        HashMap<String, HashMap<String, Object>> tmpCircularSugarMoietiesMap = new HashMap<>(2000, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                //using the clean smiles without hydrogens here because of the depiction and the hashing of molecules
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
                boolean tmpHasAnyCircularSugar = tmpSugarRemovalUtil.hasCircularSugars(tmpMolecule);
                if (tmpHasAnyCircularSugar) {
                    tmpHasCircularSugarsCounter++;
                    List<IAtomContainer> tmpRemovedMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                    if (tmpRemovedMoieties.isEmpty()) {
                        //should not happen, precaution
                        continue;
                    }
                    tmpOutput = tmpOutput.concat(tmpID + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
                    tmpOutput = tmpOutput.concat(tmpSmilesCode);
                    tmpRemovedMoieties.remove(0);
                    for (IAtomContainer tmpMoiety : tmpRemovedMoieties) {
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoiety);
                        CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpMoiety);
                        String tmpMoietySmilesCode;
                        try {
                            tmpMoietySmilesCode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException anException) {
                            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                            tmpMoietySmilesCode = "[exception]";
                        }
                        if (tmpCircularSugarMoietiesMap.containsKey(tmpMoietySmilesCode)) {
                            HashMap<String, Object> tmpInnerMap = tmpCircularSugarMoietiesMap.get(tmpMoietySmilesCode);
                            int tmpFrequency = (int)tmpInnerMap.get("FREQUENCY");
                            tmpInnerMap.put("FREQUENCY", tmpFrequency + 1);
                        } else {
                            tmpDifferentMoietiesCounter++;
                            HashMap<String, Object> tmpInnerMap = new HashMap<>(5,1);
                            tmpInnerMap.put("FREQUENCY", 1);
                            tmpInnerMap.put("FIRST_ORIGIN", tmpID);
                            tmpInnerMap.put("SMILES", tmpMoietySmilesCode);
                            tmpCircularSugarMoietiesMap.put(tmpMoietySmilesCode, tmpInnerMap);
                        }
                        tmpOutput = tmpOutput.concat(GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpMoietySmilesCode);
                    }
                    tmpCSVperMoleculeWriter.println(tmpOutput);
                    tmpCSVperMoleculeWriter.flush();
                    tmpOutput = "";
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        String[] tmpKeySet = tmpCircularSugarMoietiesMap.keySet().toArray(new String[0]);
        Arrays.sort(tmpKeySet, new Comparator<>() {
            public int compare(String aFirstKey, String aSecondKey) {
                int tmpFirstFrequency = (int)tmpCircularSugarMoietiesMap.get(aFirstKey).get("FREQUENCY");
                int tmpSecondFrequency = (int)tmpCircularSugarMoietiesMap.get(aSecondKey).get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpOutputWriter.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        System.out.println("Different detected circular sugar moieties counter: " + tmpDifferentMoietiesCounter);
        tmpOutputWriter.println("Different detected circular sugar moieties counter: " + tmpDifferentMoietiesCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Producing output CSV and images...");
        for (String tmpUniqueSmilesCode : tmpKeySet) {
            String tmpEntry = tmpUniqueSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR;
            HashMap<String, Object> tmpInnerMap = tmpCircularSugarMoietiesMap.get(tmpUniqueSmilesCode);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("FREQUENCY") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat((String)tmpInnerMap.get("FIRST_ORIGIN"));
            tmpCSVmoietyFreqWriter.println(tmpEntry);
            if ((int)tmpInnerMap.get("FREQUENCY") > 9) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpInnerMap.get("SMILES")))
                        .writeTo(tmpOutputFolderPath + File.separator + tmpInnerMap.get("FREQUENCY")
                                + "_" + tmpUniqueSmilesCode +".png");
            }
        }
        tmpOutputWriter.flush();
        tmpCSVmoietyFreqWriter.flush();
        tmpCSVperMoleculeWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVmoietyFreqWriter.close();
        tmpCSVperMoleculeWriter.close();
    }

    /**
     * This test method compiles and counts the removed linear sugar moieties the SRU returns by removing sugars from
     * all COCONUT molecules. Image files of the individual moieties with their frequencies as names are created in the
     * directory ./GlycosylationStatisticsTest_Output/coconut_stats_removed_linear_moiety_frequencies_test. Note that
     * these numbers give the total appearance of the respective moiety, not how many molecules have this moiety.
     * Also note that the linear sugar moieties returned by the SRU only consist of the detected moieties without any further
     * information/side-chains/substituents that might also be removed due to the preservation mode setting.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_removed_linear_moiety_frequencies_test.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsRemovedLinearMoietyFrequenciesTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_removed_linear_moiety_frequencies_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVmoietyFreqWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugMoietyFrequencies.csv");
        PrintWriter tmpCSVperMoleculeWriter = this.initializeOutputFile(tmpOutputFolderPath, "MoleculesWithLinSugars.csv");
        String tmpCSVperMoleculeFileHeader = "ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "originalMoleculeSMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "SugarMoietySMILES";
        tmpCSVperMoleculeWriter.println(tmpCSVperMoleculeFileHeader);
        tmpCSVperMoleculeWriter.flush();
        String tmpCSVmoietyFreqFileHeader = "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "frequency" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "firstOrigin";
        tmpCSVmoietyFreqWriter.println(tmpCSVmoietyFreqFileHeader);
        tmpCSVmoietyFreqWriter.flush();
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        String tmpOutput = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        int tmpDifferentMoietiesCounter = 0;
        HashMap<String, HashMap<String, Object>> tmpLinearSugarMoietiesMap = new HashMap<>(2000, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                //using the clean smiles without hydrogens here because of the depiction and the hashing of molecules
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                //note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
                boolean tmpHasAnyLinearSugar = tmpSugarRemovalUtil.hasLinearSugars(tmpMolecule);
                if (tmpHasAnyLinearSugar) {
                    tmpHasLinearSugarsCounter++;
                    List<IAtomContainer> tmpRemovedMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                    if (tmpRemovedMoieties.isEmpty()) {
                        //should not happen, precaution
                        continue;
                    }
                    tmpOutput = tmpOutput.concat(tmpID + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
                    tmpOutput = tmpOutput.concat(tmpSmilesCode);
                    tmpRemovedMoieties.remove(0);
                    for (IAtomContainer tmpMoiety : tmpRemovedMoieties) {
                        AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(tmpMoiety);
                        CDKHydrogenAdder.getInstance(DefaultChemObjectBuilder.getInstance()).addImplicitHydrogens(tmpMoiety);
                        String tmpMoietySmilesCode;
                        try {
                            tmpMoietySmilesCode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException anException) {
                            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                            tmpMoietySmilesCode = "[exception]";
                        }
                        if (tmpLinearSugarMoietiesMap.containsKey(tmpMoietySmilesCode)) {
                            HashMap<String, Object> tmpInnerMap = tmpLinearSugarMoietiesMap.get(tmpMoietySmilesCode);
                            int tmpFrequency = (int)tmpInnerMap.get("FREQUENCY");
                            tmpInnerMap.put("FREQUENCY", tmpFrequency + 1);
                        } else {
                            tmpDifferentMoietiesCounter++;
                            HashMap<String, Object> tmpInnerMap = new HashMap<>(5,1);
                            tmpInnerMap.put("FREQUENCY", 1);
                            tmpInnerMap.put("FIRST_ORIGIN", tmpID);
                            tmpInnerMap.put("SMILES", tmpMoietySmilesCode);
                            tmpLinearSugarMoietiesMap.put(tmpMoietySmilesCode, tmpInnerMap);
                        }
                        tmpOutput = tmpOutput.concat(GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpMoietySmilesCode);
                    }
                    tmpCSVperMoleculeWriter.println(tmpOutput);
                    tmpCSVperMoleculeWriter.flush();
                    tmpOutput = "";
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        String[] tmpKeySet = tmpLinearSugarMoietiesMap.keySet().toArray(new String[0]);
        Arrays.sort(tmpKeySet, new Comparator<>() {
            public int compare(String aFirstKey, String aSecondKey) {
                int tmpFirstFrequency = (int)tmpLinearSugarMoietiesMap.get(aFirstKey).get("FREQUENCY");
                int tmpSecondFrequency = (int)tmpLinearSugarMoietiesMap.get(aSecondKey).get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        System.out.println("Different detected linear sugar moieties counter: " + tmpDifferentMoietiesCounter);
        tmpOutputWriter.println("Different detected linear sugar moieties counter: " + tmpDifferentMoietiesCounter);
        System.out.println();
        tmpOutputWriter.println();
        System.out.println("Producing output CSV and images...");
        for (String tmpUniqueSmilesCode : tmpKeySet) {
            String tmpEntry = tmpUniqueSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR;
            HashMap<String, Object> tmpInnerMap = tmpLinearSugarMoietiesMap.get(tmpUniqueSmilesCode);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("FREQUENCY") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat((String)tmpInnerMap.get("FIRST_ORIGIN"));
            tmpCSVmoietyFreqWriter.println(tmpEntry);
            if ((int)tmpInnerMap.get("FREQUENCY") > 9) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpInnerMap.get("SMILES")))
                        .writeTo(tmpOutputFolderPath + File.separator + tmpInnerMap.get("FREQUENCY")
                                + "_" + tmpUniqueSmilesCode +".png");
            }
        }
        tmpOutputWriter.flush();
        tmpCSVmoietyFreqWriter.flush();
        tmpCSVperMoleculeWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVmoietyFreqWriter.close();
        tmpCSVperMoleculeWriter.close();
    }

    /**
     * This test method performs substructure searching and match counting in COCONUT with a set of 344 circular sugar
     * moieties reported in bacterial natural products and compiled by <a href="https://doi.org/10.1039/C4CS00426D"
     * >Elshahawi et al.</a>, additionally checking whether they are detectable by the SRU. The structure are first
     * grouped by stereoisomers. Then, they were separated in two datasets, depending on whether they are detectable
     * with the SRU default options as sugar moieties or not. Finally, a substructure search for every stereoisomer
     * group of both datasets is performed in COCONUT, and the matches counted. Note that the number of substructure
     * matches detected here would not precisely correspond to the number of cases where the SRU would remove the
     * glycosidic moiety as the SRU does not only detect the presence of a substructure but also detects whether it is
     * terminal and/or an isolated cycle.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_review_data_sugars_test. Additionally, CSV files
     * and image files of the most frequent SRU-positive and -negative moieties are created, respectively.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsReviewDataSugarsAppearanceTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_review_data_sugars_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVSRUPositiveWriter = this.initializeOutputFile(tmpOutputFolderPath, "SRUSugars.csv");
        PrintWriter tmpCSVSRUNegativeWriter = this.initializeOutputFile(tmpOutputFolderPath, "NonSRUSugars.csv");
        tmpCSVSRUPositiveWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpCSVSRUNegativeWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = null;
        try {
            tmpSDFile = new File(tmpClassLoader.getResource("review_glycosylated_NPs_bacteria_data.sdf").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("SDF could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpSDFile.getAbsolutePath());
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), true);
        HashMap<String, HashMap<String, Object>> tmpSRUPositiveSugarPatterns = new HashMap<>(344, 1);
        HashMap<String, HashMap<String, Object>> tmpSRUNegativeSugarPatterns = new HashMap<>(344, 1);
        String tmpReviewSugarID;
        int tmpReviewSugarsCounter = 0;
        int tmpReviewDataExceptionsCounter = 0;
        while (tmpReader.hasNext()) {
            tmpReviewSugarID = "[unidentified]";
            try {
                IAtomContainer tmpReviewSugar = tmpReader.next();
                tmpReviewSugarsCounter++;
                tmpReviewSugarID = tmpReviewSugar.getProperty("Name");
                HashMap<String, Object> tmpMap = new HashMap<>(4, 1);
                tmpMap.put("PATTERN", DfPattern.findSubstructure(tmpReviewSugar));
                tmpMap.put("ID", tmpReviewSugarID);
                tmpMap.put("FREQUENCY", 0);
                String tmpSmilesCode = "[generation_failed]";
                try {
                    tmpSmilesCode = tmpSmiGen.create(tmpReviewSugar);
                } catch (CDKException anException) {
                    GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpReviewSugarID, anException);
                    tmpReviewDataExceptionsCounter++;
                    continue;
                }
                tmpMap.put("SMILES", tmpSmilesCode);
                boolean tmpHasSugars = tmpSugarRemovalUtil.hasCircularSugars(tmpReviewSugar);
                if (tmpHasSugars) {
                    if (tmpSRUPositiveSugarPatterns.containsKey(tmpSmilesCode)) {
                        HashMap<String, Object> tmpInnerMap = tmpSRUPositiveSugarPatterns.get(tmpSmilesCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUPositiveSugarPatterns.put(tmpSmilesCode, tmpMap);
                    }
                } else {
                    if (tmpSRUNegativeSugarPatterns.containsKey(tmpSmilesCode)) {
                        HashMap<String, Object> tmpInnerMap = tmpSRUNegativeSugarPatterns.get(tmpSmilesCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUNegativeSugarPatterns.put(tmpSmilesCode, tmpMap);
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpReviewSugarID, anException);
                tmpReviewDataExceptionsCounter++;
            }
        }
        System.out.println("Parsing of patterns done, iterating COCONUT...");
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUPositiveSugarPatterns.values()) {
                    DfPattern tmpPattern = (DfPattern) tmpReviewSugarMap.get("PATTERN");
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpReviewSugarMap.put("FREQUENCY", ((int)tmpReviewSugarMap.get("FREQUENCY") + 1));
                    }
                }
                for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUNegativeSugarPatterns.values()) {
                    DfPattern tmpPattern = (DfPattern) tmpReviewSugarMap.get("PATTERN");
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpReviewSugarMap.put("FREQUENCY", ((int)tmpReviewSugarMap.get("FREQUENCY") + 1));
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions in the review data counter: " + tmpReviewDataExceptionsCounter);
        tmpOutputWriter.println("Exceptions in the review data counter: " + tmpReviewDataExceptionsCounter);
        System.out.println("Molecules in the review data counter: " + tmpReviewSugarsCounter);
        tmpOutputWriter.println("Molecules in the review data counter: " + tmpReviewSugarsCounter);
        System.out.println("Exceptions in COCONUT counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions in COCONUT counter: " + tmpExceptionsCounter);
        System.out.println("Molecules in COCONUT counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules in COCONUT counter: " + tmpMoleculesCounter);
        int tmpNumberOfDistinctPatterns = tmpSRUPositiveSugarPatterns.size() + tmpSRUNegativeSugarPatterns.size();
        System.out.println("How many distinct moieties remained after grouping of stereo-isomers: " + tmpNumberOfDistinctPatterns);
        tmpOutputWriter.println("How many distinct moieties remained after grouping of stereo-isomers: " + tmpNumberOfDistinctPatterns);
        System.out.println(tmpSRUPositiveSugarPatterns.size() + " of these were detectable by the SRU.");
        tmpOutputWriter.println(tmpSRUPositiveSugarPatterns.size() + " of these were detectable by the SRU.");
        System.out.println(tmpSRUNegativeSugarPatterns.size() + " of these were NOT detectable by the SRU.");
        tmpOutputWriter.println(tmpSRUNegativeSugarPatterns.size() + " of these were NOT detectable by the SRU.");
        HashMap<String, Object>[] tmpSRUPositivePatternsMapArray = tmpSRUPositiveSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUPositivePatternsMapArray, new Comparator<>() {
            public int compare(HashMap aFirstMap, HashMap aSecondMap) {
                int tmpFirstFrequency = (int)aFirstMap.get("FREQUENCY");
                int tmpSecondFrequency = (int)aSecondMap.get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        HashMap<String, Object>[] tmpSRUNegativePatternsMapArray = tmpSRUNegativeSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUNegativePatternsMapArray, new Comparator<>() {
            public int compare(HashMap aFirstMap, HashMap aSecondMap) {
                int tmpFirstFrequency = (int)aFirstMap.get("FREQUENCY");
                int tmpSecondFrequency = (int)aSecondMap.get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        File tmpSRUPositiveOutputFolder = new File (tmpOutputFolderPath + File.separator + "SRUPositive" + File.separator);
        if (!tmpSRUPositiveOutputFolder.exists()) {
            tmpSRUPositiveOutputFolder.mkdirs();
        }
        for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUPositivePatternsMapArray) {
            String tmpName = (String) tmpReviewSugarMap.get("ID");
            String tmpReviewSugarSmilesCode = (String) tmpReviewSugarMap.get("SMILES");
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUPositiveWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpReviewSugarSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && !(tmpReviewSugarSmilesCode.equals("[generation_failed]"))) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles(tmpReviewSugarSmilesCode))
                        .writeTo(tmpOutputFolderPath + File.separator + "SRUPositive" + File.separator + tmpFrequency
                                + "_" + tmpName +".png");
            }
        }
        File tmpSRUNegativeOutputFolder = new File (tmpOutputFolderPath + File.separator + "SRUNegative" + File.separator);
        if (!tmpSRUNegativeOutputFolder.exists()) {
            tmpSRUNegativeOutputFolder.mkdirs();
        }
        for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUNegativePatternsMapArray) {
            String tmpName = (String) tmpReviewSugarMap.get("ID");
            String tmpReviewSugarSmilesCode = (String) tmpReviewSugarMap.get("SMILES");
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUNegativeWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpReviewSugarSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && !(tmpReviewSugarSmilesCode.equals("[generation_failed]"))) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles(tmpReviewSugarSmilesCode))
                        .writeTo(tmpOutputFolderPath + File.separator + "SRUNegative" + File.separator + tmpFrequency
                                + "_" + tmpName +".png");
            }
        }
        tmpOutputWriter.flush();
        tmpCSVSRUPositiveWriter.flush();
        tmpCSVSRUNegativeWriter.flush();
        tmpReader.close();
        tmpOutputWriter.close();
        tmpCSVSRUPositiveWriter.close();
        tmpCSVSRUNegativeWriter.close();
    }

    /**
     * This test method performs substructure searching and match counting in COCONUT with a set of 344 circular sugar
     * moieties reported in bacterial natural products and compiled by <a href="https://doi.org/10.1039/C4CS00426D"
     * >Elshahawi et al.</a>, additionally checking whether they are detectable by the SRU. In contrast to the other
     * test method doing the same analysis, here, the exocyclic oxygen ratio threshold setting of the SRU is lowered to 0.3
     * instead of using the default value of 0.5. The structure are first
     * grouped by stereoisomers. Then, they were separated in two datasets, depending on whether they are detectable
     * with the SRU default options as sugar moieties or not. Finally, a substructure search for every stereoisomer
     * group of both datasets is performed in COCONUT, and the matches counted. Note that the number of substructure
     * matches detected here would not precisely correspond to the number of cases where the SRU would remove the
     * glycosidic moiety as the SRU does not only detect the presence of a substructure but also detects whether it is
     * terminal and/or an isolated cycle.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/coconut_stats_review_data_sugars_nondefault_test. Additionally, CSV files
     * and image files of the most frequent SRU-positive and -negative moieties are created, respectively.
     * Test is ignored, if no connection to MongoDB can be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutStatsReviewDataSugarsAppearanceNonDefaultTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        final String tmpSpecificOutputFolderName = "coconut_stats_review_data_sugars_nondefault_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpCSVSRUPositiveWriter = this.initializeOutputFile(tmpOutputFolderPath, "SRUSugars.csv");
        PrintWriter tmpCSVSRUNegativeWriter = this.initializeOutputFile(tmpOutputFolderPath, "NonSRUSugars.csv");
        tmpCSVSRUPositiveWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpCSVSRUNegativeWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = null;
        try {
            tmpSDFile = new File(tmpClassLoader.getResource("review_glycosylated_NPs_bacteria_data.sdf").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("SDF could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpSDFile.getAbsolutePath());
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //change setting to detect more of the sugar-like moieties
        tmpSugarRemovalUtil.setExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting(0.3);
        System.out.println("SRU exocyclic oxygen to atoms in ring ratio threshold set to "
                + tmpSugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        tmpOutputWriter.println("SRU exocyclic oxygen to atoms in ring ratio threshold set to "
                + tmpSugarRemovalUtil.getExocyclicOxygenAtomsToAtomsInRingRatioThresholdSetting());
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), true);
        HashMap<String, HashMap<String, Object>> tmpSRUPositiveSugarPatterns = new HashMap<>(344, 1);
        HashMap<String, HashMap<String, Object>> tmpSRUNegativeSugarPatterns = new HashMap<>(344, 1);
        String tmpReviewSugarID;
        int tmpReviewSugarsCounter = 0;
        int tmpReviewDataExceptionsCounter = 0;
        while (tmpReader.hasNext()) {
            tmpReviewSugarID = "[unidentified]";
            try {
                IAtomContainer tmpReviewSugar = tmpReader.next();
                tmpReviewSugarsCounter++;
                tmpReviewSugarID = tmpReviewSugar.getProperty("Name");
                HashMap<String, Object> tmpMap = new HashMap<>(4, 1);
                tmpMap.put("PATTERN", DfPattern.findSubstructure(tmpReviewSugar));
                tmpMap.put("ID", tmpReviewSugarID);
                tmpMap.put("FREQUENCY", 0);
                String tmpSmilesCode = "[generation_failed]";
                try {
                    tmpSmilesCode = tmpSmiGen.create(tmpReviewSugar);
                } catch (CDKException anException) {
                    GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpReviewSugarID, anException);
                    tmpReviewDataExceptionsCounter++;
                    continue;
                }
                tmpMap.put("SMILES", tmpSmilesCode);
                boolean tmpHasSugars = tmpSugarRemovalUtil.hasCircularSugars(tmpReviewSugar);
                if (tmpHasSugars) {
                    if (tmpSRUPositiveSugarPatterns.containsKey(tmpSmilesCode)) {
                        HashMap<String, Object> tmpInnerMap = tmpSRUPositiveSugarPatterns.get(tmpSmilesCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUPositiveSugarPatterns.put(tmpSmilesCode, tmpMap);
                    }
                } else {
                    if (tmpSRUNegativeSugarPatterns.containsKey(tmpSmilesCode)) {
                        HashMap<String, Object> tmpInnerMap = tmpSRUNegativeSugarPatterns.get(tmpSmilesCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUNegativeSugarPatterns.put(tmpSmilesCode, tmpMap);
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpReviewSugarID, anException);
                tmpReviewDataExceptionsCounter++;
            }
        }
        System.out.println("Parsing of patterns done, iterating COCONUT...");
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUPositiveSugarPatterns.values()) {
                    DfPattern tmpPattern = (DfPattern) tmpReviewSugarMap.get("PATTERN");
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpReviewSugarMap.put("FREQUENCY", ((int)tmpReviewSugarMap.get("FREQUENCY") + 1));
                    }
                }
                for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUNegativeSugarPatterns.values()) {
                    DfPattern tmpPattern = (DfPattern) tmpReviewSugarMap.get("PATTERN");
                    if (tmpPattern.matches(tmpMolecule)) {
                        tmpReviewSugarMap.put("FREQUENCY", ((int)tmpReviewSugarMap.get("FREQUENCY") + 1));
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions in the review data counter: " + tmpReviewDataExceptionsCounter);
        tmpOutputWriter.println("Exceptions in the review data counter: " + tmpReviewDataExceptionsCounter);
        System.out.println("Molecules in the review data counter: " + tmpReviewSugarsCounter);
        tmpOutputWriter.println("Molecules in the review data counter: " + tmpReviewSugarsCounter);
        System.out.println("Exceptions in COCONUT counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions in COCONUT counter: " + tmpExceptionsCounter);
        System.out.println("Molecules in COCONUT counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules in COCONUT counter: " + tmpMoleculesCounter);
        int tmpNumberOfDistinctPatterns = tmpSRUPositiveSugarPatterns.size() + tmpSRUNegativeSugarPatterns.size();
        System.out.println("How many distinct moieties remained after grouping of stereo-isomers: " + tmpNumberOfDistinctPatterns);
        tmpOutputWriter.println("How many distinct moieties remained after grouping of stereo-isomers: " + tmpNumberOfDistinctPatterns);
        System.out.println(tmpSRUPositiveSugarPatterns.size() + " of these were detectable by the SRU.");
        tmpOutputWriter.println(tmpSRUPositiveSugarPatterns.size() + " of these were detectable by the SRU.");
        System.out.println(tmpSRUNegativeSugarPatterns.size() + " of these were NOT detectable by the SRU.");
        tmpOutputWriter.println(tmpSRUNegativeSugarPatterns.size() + " of these were NOT detectable by the SRU.");
        HashMap<String, Object>[] tmpSRUPositivePatternsMapArray = tmpSRUPositiveSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUPositivePatternsMapArray, new Comparator<>() {
            public int compare(HashMap aFirstMap, HashMap aSecondMap) {
                int tmpFirstFrequency = (int)aFirstMap.get("FREQUENCY");
                int tmpSecondFrequency = (int)aSecondMap.get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        HashMap<String, Object>[] tmpSRUNegativePatternsMapArray = tmpSRUNegativeSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUNegativePatternsMapArray, new Comparator<>() {
            public int compare(HashMap aFirstMap, HashMap aSecondMap) {
                int tmpFirstFrequency = (int)aFirstMap.get("FREQUENCY");
                int tmpSecondFrequency = (int)aSecondMap.get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        File tmpSRUPositiveOutputFolder = new File (tmpOutputFolderPath + File.separator + "SRUPositive" + File.separator);
        if (!tmpSRUPositiveOutputFolder.exists()) {
            tmpSRUPositiveOutputFolder.mkdirs();
        }
        for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUPositivePatternsMapArray) {
            String tmpName = (String) tmpReviewSugarMap.get("ID");
            String tmpReviewSugarSmilesCode = (String) tmpReviewSugarMap.get("SMILES");
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUPositiveWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpReviewSugarSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && !(tmpReviewSugarSmilesCode.equals("[generation_failed]"))) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles(tmpReviewSugarSmilesCode))
                        .writeTo(tmpOutputFolderPath + File.separator + "SRUPositive" + File.separator + tmpFrequency
                                + "_" + tmpName +".png");
            }
        }
        File tmpSRUNegativeOutputFolder = new File (tmpOutputFolderPath + File.separator + "SRUNegative" + File.separator);
        if (!tmpSRUNegativeOutputFolder.exists()) {
            tmpSRUNegativeOutputFolder.mkdirs();
        }
        for (HashMap<String, Object> tmpReviewSugarMap : tmpSRUNegativePatternsMapArray) {
            String tmpName = (String) tmpReviewSugarMap.get("ID");
            String tmpReviewSugarSmilesCode = (String) tmpReviewSugarMap.get("SMILES");
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUNegativeWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpReviewSugarSmilesCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && !(tmpReviewSugarSmilesCode.equals("[generation_failed]"))) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles(tmpReviewSugarSmilesCode))
                        .writeTo(tmpOutputFolderPath + File.separator + "SRUNegative" + File.separator + tmpFrequency
                                + "_" + tmpName +".png");
            }
        }
        tmpOutputWriter.flush();
        tmpCSVSRUPositiveWriter.flush();
        tmpCSVSRUNegativeWriter.flush();
        tmpReader.close();
        tmpOutputWriter.close();
        tmpCSVSRUPositiveWriter.close();
        tmpCSVSRUNegativeWriter.close();
    }

    /**
     * This test method does some basic glycosylation statistics of COCONUT. The important thing here is that it uses
     * the COCONUT SDF instead of connecting to a MongoDB COCONUT instance. the method calculates: How many molecules have sugars,
     * circular sugars, linear sugars, both, and how many molecules are basically sugars. All statistics are printed to
     * console and also compiled in an output file created in the directory ./GlycosylationStatisticsTest_Output/coconut_sdf_test/.
     * Test is ignored, if the COCONUT SDF cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void coconutSdfTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = null;
        try {
            tmpSDFile = new File(tmpClassLoader.getResource(GlycosylationStatisticsTest.SDF_NAME).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("COCONUT SDF could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpSDFile.getAbsolutePath());
        final String tmpSpecificOutputFolderName = "coconut_sdf_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), true);
        String tmpID;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsCNPs = new ArrayList<>(50000);
        int tmpHasNoSugarsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList<>(2000);
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarCNPs = new ArrayList<>(2000);
        while (tmpReader.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpMolecule = tmpReader.next();
                tmpMoleculesCounter++;
                tmpID = tmpMolecule.getProperty(GlycosylationStatisticsTest.ID_KEY);
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not detected/removed/counted!
                // note also: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                boolean tmpHasAnyTypeOfSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyCircularSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyLinearSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpHasAnyTypeOfSugarsCounter++;
                    tmpHasAnyTypeOfSugarsCNPs.add(tmpID);
                    if (tmpHasAnyCircularSugar) {
                        tmpHasCircularSugarsCounter++;
                        tmpHasCircularSugarsCNPs.add(tmpID);
                    }
                    if (tmpHasAnyLinearSugar) {
                        tmpHasLinearSugarsCounter++;
                        tmpHasLinearSugarsCNPs.add(tmpID);
                    }
                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsCNPs.add(tmpID);
                    }
                    if (tmpDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                        tmpBasicallyASugarCNPs.add(tmpID);
                    }
                } else {
                    tmpHasNoSugarsCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        tmpOutputWriter.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        double tmpPercentage = ((double) tmpHasAnyTypeOfSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain sugars.");
        System.out.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        tmpOutputWriter.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpOutputWriter.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        tmpOutputWriter.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Sugar-containing molecules: " + tmpHasAnyTypeOfSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Circular-sugar-containing molecules: " + tmpHasCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Linear-sugar-containing molecules: " + tmpHasLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules containing both circular and linear sugars: " + tmpHasCircularAndLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a sugar: " + tmpBasicallyASugarCNPs);
        tmpOutputWriter.flush();
        tmpReader.close();
        tmpOutputWriter.close();
        Assert.assertEquals(tmpMoleculesCounter, tmpHasNoSugarsCounter + tmpHasAnyTypeOfSugarsCounter);
    }

    /**
     * This test method checks whether there are molecules in the COCONUT MongoDB instance that are missing in the COCONUT
     * SDF. All statistics are printed to console. Test is ignored, if the COCONUT SDF cannot be found or no connection
     * to MongoDB be made.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void checkCoconutSdfTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = null;
        try {
            tmpSDFile = new File(tmpClassLoader.getResource(GlycosylationStatisticsTest.SDF_NAME).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("COCONUT SDF could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpSDFile.getAbsolutePath());
        //Note: Skip is set to false!!
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), false);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        HashSet<String> tmpIDsInSDFset = new HashSet<>(450000);
        String tmpID;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpReader.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpMolecule = tmpReader.next();
                tmpMoleculesCounter++;
                tmpID = tmpMolecule.getProperty(GlycosylationStatisticsTest.ID_KEY);
                tmpIDsInSDFset.add(tmpID);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println("Done with the SDF.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Now analyzing the MongoDB database and searching for missing IDs...");
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        Document tmpCurrentDoc;
        String tmpSmilesCode;
        tmpMoleculesCounter = 0;
        tmpExceptionsCounter = 0;
        int tmpMissingMoleculesCounter = 0;
        System.out.println("Following CNP IDs are missing in the SDF:");
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                if (!tmpIDsInSDFset.contains(tmpID)) {
                    System.out.println(tmpID);
                    tmpMissingMoleculesCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Count of molecules missing in the SDF: " + tmpMissingMoleculesCounter);
        tmpCursor.close();
    }
    //</editor-fold>

    //<editor-fold desc="ZINC">
    /**
     * This test method calculates general glycosylation statistics for a curated set of synthetic molecules taken from the
     * ZINC "for-sale" subset. 475958 molecules are contained in the curated dataset. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/zinc_for-sale_stats_basics_test.
     * Test is ignored, if specified dataset cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincSyntheticsForSaleStatsBasicsTest() throws Exception {
        this.doBasicGlycoStatsOnSMILESfile("ZINC_for-sale_without_biogenics_and_COCONUT.txt",
                "zinc_for-sale_stats_basics_test");
    }

    /**
     * This test method calculates general glycosylation statistics for a curated set of synthetic molecules taken from the
     * ZINC "in-vitro" subset. For these molecules, activities in in vitro assays have been reported. 65620 molecules
     * are contained in the curated dataset. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/zinc_in-vitro_stats_basics_test.
     * Test is ignored, if specified dataset cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincSyntheticsInVitroStatsBasicsTest() throws Exception {
        this.doBasicGlycoStatsOnSMILESfile("ZINC_in-vitro_without_biogenics_and_COCONUT.txt",
                "zinc_in-vitro_stats_basics_test");
    }

    /**
     * This test method calculates general glycosylation statistics for a curated set of molecules taken from the
     * ZINC "in-vitro" subset. For these molecules, activities in in vitro assays have been reported. 306347 molecules
     * are contained in the curated dataset. The only curation step here was to keep only one molecule of identified
     * stereo-isomers, respectively. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/zinc_in-vitro_complete_stats_basics_test.
     * Test is ignored, if specified dataset cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincCompleteInVitroSubsetStatsBasicsTest() throws Exception {
        this.doBasicGlycoStatsOnSMILESfile("ZINC_in-vitro_flattened.txt",
                "zinc_in-vitro_complete_stats_basics_test");
    }

    /**
     * Curation test method for ZINC "synthetics for-sale" dataset. The structures in the given SMILES file are iterated
     * and matched to all COCONUT and ZINC "biogenic" subset molecules based on unique SMILES codes. All matching
     * molecules are filtered to obtain a dataset of "synthetic" molecules. Additionally, only one structure of every
     * stereo-isomer group in the original dataset is retained. The curated set is written to file as SMILES codes in
     * the specified output folder ./GlycosylationStatisticsTest_Output/zinc_for-sale_curation_test.
     * To analyse the curated dataset, it needs to be put into the test resource directory!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincForSaleDatasetCurationTest() throws Exception {
        this.filterZINCBiogenicAndCOCONUTMoleculesAndGroupStereoIsomersOfSMILESfile(
                "ZINC_for-sale_picked_subset_2020_Okt_19.txt",
                500000,
                "ZINC_for-sale_without_biogenics_and_COCONUT.txt",
                "zinc_for-sale_curation_test");
    }

    /**
     * Curation test method for ZINC "synthetics in-vitro" dataset. The structures in the given SMILES file are iterated
     * and matched to all COCONUT and ZINC "biogenic" subset molecules based on unique SMILES codes. All matching
     * molecules are filtered to obtain a dataset of "synthetic" molecules. Additionally, only one structure of every
     * stereo-isomer group in the original dataset is retained. The curated set is written to file as SMILES codes in
     * the specified output folder ./GlycosylationStatisticsTest_Output/zinc_in-vitro_curation_test.
     * To analyse the curated dataset, it needs to be put into the test resource directory!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincInVitroDatasetCurationTest() throws Exception {
        this.filterZINCBiogenicAndCOCONUTMoleculesAndGroupStereoIsomersOfSMILESfile(
                "ZINC_in-vitro_subset_2020_Okt_30.txt",
                306347,
                "ZINC_in-vitro_without_biogenics_and_COCONUT.txt",
                "zinc_in-vitro_curation_test");
    }

    /**
     * Curation test method for ZINC "in-vitro" dataset. From the structures in the given SMILES file, only one structure
     * of every stereo-isomer group is retained. The curated set is written to file as SMILES codes in
     * the specified output folder ./GlycosylationStatisticsTest_Output/zinc_in-vitro_complete_curation_test.
     * To analyse the curated dataset, it needs to be put into the test resource directory!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void zincInVitroCompleteDatasetCurationTest() throws Exception {
        this.groupStereoisomersOfSMILESfile("ZINC_in-vitro_subset_2020_Okt_30.txt", 306347,
                "ZINC_in-vitro_flattened.txt", "zinc_in-vitro_complete_curation_test");
    }
    //</editor-fold>
    //<editor-fold desc="ChEMBL">
    /**
     * Curation test method for ChEMBL dataset. From the structures in the given SDF, only one structure
     * of every stereo-isomer group is retained. The curated set is written to file as SMILES codes in
     * the specified output folder ./GlycosylationStatisticsTest_Output/chembl_curation_test.
     * To analyse the curated dataset, it needs to be put into the test resource directory!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void chemblCurationTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger("chembl_curation_test");
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //Loading analysed dataset
        System.out.println("Loading and processing the given dataset now...");
        tmpOutputWriter.println("Loading and processing the given dataset now...");
        File tmpOriginalDataSetSmilesFile = null;
        try {
            tmpOriginalDataSetSmilesFile = new File(tmpClassLoader.getResource("chembl_28.sdf").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("Original data set file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Data set found at: " + tmpOriginalDataSetSmilesFile.getAbsolutePath());
        FileReader tmpOriginalDataSetSmilesFileReader = new FileReader(tmpOriginalDataSetSmilesFile);
        BufferedReader tmpOriginalDataSetSmilesBufferedReader = new BufferedReader(tmpOriginalDataSetSmilesFileReader);
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(tmpOriginalDataSetSmilesBufferedReader, DefaultChemObjectBuilder.getInstance(), true);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique); //Unique does not encode stereochemistry! Absolute would do
        //key: unique SMILES string (does not encode stereochemistry!)
        //object: HashMap of ID of the first stereo isomer encountered and number of stereo-isomers encountered
        HashMap<String, HashMap<String, Object>> tmpSmilesMap = new HashMap<>((int)(2066376*1.1), 1.0f);
        String tmpID = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpSDFReader.hasNext()) {
            try {
                IAtomContainer tmpMoleculeFromFile = tmpSDFReader.next();
                if (Objects.isNull(tmpMoleculeFromFile)) {
                    break;
                }
                tmpMoleculesCounter++;
                if ((tmpMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpMoleculesCounter + " molecules were processed already...");
                }
                tmpID = tmpMoleculeFromFile.getProperty("chembl_id");
                String tmpCDKSMILESCode = tmpSmiGen.create(tmpMoleculeFromFile);
                if (tmpSmilesMap.containsKey(tmpCDKSMILESCode)) {
                    HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpCDKSMILESCode);
                    int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
                    tmpInnerMap.put("FREQUENCY", (tmpStereoIsomerNumber + 1));
                    continue;
                } else {
                    HashMap<String, Object> tmpInnerMap = new HashMap<>(3, 1.0f);
                    tmpInnerMap.put("ID", tmpID);
                    tmpInnerMap.put("FREQUENCY", 1);
                    tmpSmilesMap.put(tmpCDKSMILESCode, tmpInnerMap);
                    continue;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpOriginalDataSetSmilesFileReader.close();
        tmpOriginalDataSetSmilesBufferedReader.close();
        System.out.println("Processing of the given dataset done.");
        tmpOutputWriter.println("Processing of the given dataset done.");
        System.out.println(tmpMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        //Writing to file
        System.out.println("Writing the curated data set to file now...");
        PrintWriter tmpMoleculesWriter = this.initializeOutputFile(tmpOutputFolderPath, "chembl_28_curated.txt");
        int tmpCounter = 0;
        tmpID = "";
        for (String tmpSmilesCode : tmpSmilesMap.keySet()) {
            if ((tmpCounter % 10000) == 0) {
                System.out.println(tmpCounter + " lines were processed already...");
            }
            HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpSmilesCode);
            int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
            tmpID = (String) tmpInnerMap.get("ID");
            tmpMoleculesWriter.println(tmpSmilesCode + " " + tmpID + " " + tmpStereoIsomerNumber);
            tmpCounter++;
        }
        System.out.println("Done, shutting down.");
        tmpMoleculesWriter.flush();
        tmpMoleculesWriter.close();
        tmpOutputWriter.flush();
        tmpOutputWriter.close();
    }

    /**
     * This test method calculates general glycosylation statistics for a curated set of molecules taken from the
     * ChEMBL database. 1982219 molecules are contained in the curated dataset. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/chembl_stats_test.
     * Test is ignored, if specified dataset cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void chemblStatisticsTest() throws Exception {
        this.doBasicGlycoStatsOnSMILESfile("chembl_28_curated.txt", "chembl_stats_test");
    }
    //</editor-fold>
    //
    //<editor-fold desc="DrugBank">
    /**
     * Curation test method for DrugBank dataset. From the structures in the given SDF, only one structure
     * of every stereo-isomer group is retained. The curated set is written to file as SMILES codes in
     * the specified output folder ./GlycosylationStatisticsTest_Output/chembl_curation_test.
     * To analyse the curated dataset, it needs to be put into the test resource directory!
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void drugBankCurationTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger("drugbank_curation_test");
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //Loading analysed dataset
        System.out.println("Loading and processing the given dataset now...");
        tmpOutputWriter.println("Loading and processing the given dataset now...");
        File tmpOriginalDataSetSmilesFile = null;
        try {
            tmpOriginalDataSetSmilesFile = new File(tmpClassLoader.getResource("drugbank_all_structures_2021_03_11.sdf").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("Original data set file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Data set found at: " + tmpOriginalDataSetSmilesFile.getAbsolutePath());
        FileReader tmpOriginalDataSetSmilesFileReader = new FileReader(tmpOriginalDataSetSmilesFile);
        BufferedReader tmpOriginalDataSetSmilesBufferedReader = new BufferedReader(tmpOriginalDataSetSmilesFileReader);
        IteratingSDFReader tmpSDFReader = new IteratingSDFReader(tmpOriginalDataSetSmilesBufferedReader, DefaultChemObjectBuilder.getInstance(), true);
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique); //Unique does not encode stereochemistry! Absolute would do
        //key: unique SMILES string (does not encode stereochemistry!)
        //object: HashMap of ID of the first stereo isomer encountered and number of stereo-isomers encountered
        HashMap<String, HashMap<String, Object>> tmpSmilesMap = new HashMap<>((int)(11172*1.1), 1.0f);
        String tmpID = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (tmpSDFReader.hasNext()) {
            try {
                IAtomContainer tmpMoleculeFromFile = tmpSDFReader.next();
                if (Objects.isNull(tmpMoleculeFromFile)) {
                    break;
                }
                tmpMoleculesCounter++;
                if ((tmpMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpMoleculesCounter + " molecules were processed already...");
                }
                tmpID = tmpMoleculeFromFile.getProperty("DATABASE_ID");
                String tmpCDKSMILESCode = tmpSmiGen.create(tmpMoleculeFromFile);
                if (tmpSmilesMap.containsKey(tmpCDKSMILESCode)) {
                    HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpCDKSMILESCode);
                    int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
                    tmpInnerMap.put("FREQUENCY", (tmpStereoIsomerNumber + 1));
                    continue;
                } else {
                    HashMap<String, Object> tmpInnerMap = new HashMap<>(3, 1.0f);
                    tmpInnerMap.put("ID", tmpID);
                    tmpInnerMap.put("FREQUENCY", 1);
                    tmpSmilesMap.put(tmpCDKSMILESCode, tmpInnerMap);
                    continue;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpOriginalDataSetSmilesFileReader.close();
        tmpOriginalDataSetSmilesBufferedReader.close();
        System.out.println("Processing of the given dataset done.");
        tmpOutputWriter.println("Processing of the given dataset done.");
        System.out.println(tmpMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        //Writing to file
        System.out.println("Writing the curated data set to file now...");
        PrintWriter tmpMoleculesWriter = this.initializeOutputFile(tmpOutputFolderPath, "drugbank_curated.txt");
        int tmpCounter = 0;
        tmpID = "";
        for (String tmpSmilesCode : tmpSmilesMap.keySet()) {
            if ((tmpCounter % 10000) == 0) {
                System.out.println(tmpCounter + " lines were processed already...");
            }
            HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpSmilesCode);
            int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
            tmpID = (String) tmpInnerMap.get("ID");
            tmpMoleculesWriter.println(tmpSmilesCode + " " + tmpID + " " + tmpStereoIsomerNumber);
            tmpCounter++;
        }
        System.out.println("Done, shutting down.");
        tmpMoleculesWriter.flush();
        tmpMoleculesWriter.close();
        tmpOutputWriter.flush();
        tmpOutputWriter.close();
    }

    /**
     * This test method calculates general glycosylation statistics for a curated set of molecules taken from the
     * DrugBank database. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the directory
     * ./GlycosylationStatisticsTest_Output/drugbank_statistics_test.
     * Test is ignored, if specified dataset cannot be found.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void drugbankStatisticsTest() throws Exception {
        this.doBasicGlycoStatsOnSMILESfile("drugbank_curated.txt", "drugbank_statistics_test");
    }
    //</editor-fold>
    //</editor-fold>
    /**
     * Generates depictions to illustrate the deglycosylation of one example molecule. Images files are created of
     * the original molecule, the original molecule with highlighted linear sugar candidates (all together and separately)
     * and the remaining aglycon after deglycosylation. All image files are created in the directory
     * ./GlycosylationStatisticsTest_Output/single_depiction_test.
     *
     * @throws Exception if anything goes wrong
     */
    @Test
    public void singleSugarDetectionDepictionTest() throws Exception {
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger("single_depiction_test");
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        IAtomContainer tmpOriginalMolecule;
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
        tmpOriginalMolecule = tmpSmiPar.parseSmiles(
                "O=C[C@H](O)[C@@H](O)[C@H](O)[C@H](O)COS(=O)(=O)O");
        tmpDepictionGenerator.withSize(2000, 2000)
                .withFillToFit()
                .depict(tmpOriginalMolecule)
                .writeTo(tmpOutputFolderPath + File.separator + "Test_original_molecule.png");
        List<IAtomContainer> tmpCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpOriginalMolecule);
        List<IAtomContainer> tmpToHighlight = new ArrayList<>(tmpCandidates.size());
        for (int i = 0; i < tmpCandidates.size(); i++) {
            IAtomContainer tmpCandidate = tmpCandidates.get(i);
            tmpToHighlight.add(tmpCandidate);
        }
        tmpDepictionGenerator.withHighlight(tmpToHighlight, Color.BLUE)
                .withSize(2000, 2000)
                .withFillToFit()
                .depict(tmpOriginalMolecule)
                .writeTo(tmpOutputFolderPath + File.separator + "Test_all_candidates" + ".png");
        for (int i = 0; i < tmpCandidates.size(); i++) {
            IAtomContainer tmpCandidate = tmpCandidates.get(i);
            List<IAtomContainer> tmpCandidateList = new ArrayList<>(1);
            tmpCandidateList.add(tmpCandidate);
            tmpDepictionGenerator.withHighlight(tmpCandidateList, Color.BLUE)
                    .withSize(2000, 2000)
                    .withFillToFit()
                    .depict(tmpOriginalMolecule)
                    .writeTo(tmpOutputFolderPath + File.separator + "Test_candidates_separately_" + i + ".png");
        }
        IAtomContainer tmpAglycon = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpOriginalMolecule, true);
        tmpDepictionGenerator.withSize(2000, 2000)
                .withFillToFit()
                .depict(tmpAglycon)
                .writeTo(tmpOutputFolderPath + File.separator + "Test_deglycosylated_molecule.png");
    }
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * Returns a MongoCursor object for iterating the COCONUT database, using the credentials specified in the
     * respective private static final constants. Throws MongoTimeoutException if no connection can be made.
     */
    private MongoCursor<Document> getCOCONUTMongoCursorForIteration() throws MongoTimeoutException {
        MongoClientSettings.Builder tmpBuilder = MongoClientSettings.builder();
        ServerAddress tmpAddress = new ServerAddress(GlycosylationStatisticsTest.HOST, GlycosylationStatisticsTest.PORT);
        tmpBuilder.applyToClusterSettings(builder -> builder.hosts(Collections.singletonList(tmpAddress)));
        MongoClientSettings tmpSettings = tmpBuilder.build();
        MongoClient tmpMongoClient = MongoClients.create(tmpSettings);
        MongoDatabase tmpDatabase = tmpMongoClient.getDatabase(GlycosylationStatisticsTest.DATABASE_NAME);
        MongoCollection<Document> tmpCollection = tmpDatabase.getCollection(GlycosylationStatisticsTest.COLLECTION_NAME);
        MongoCursor<Document> tmpCursor;
        tmpCursor = tmpCollection.find().iterator();
        System.out.println("Connection to MongoDB successful.");
        System.out.println("Collection " + tmpCollection.getNamespace().getCollectionName() + " in database " + tmpDatabase.getName() + " is loaded.");
        return tmpCursor;
    }

    /**
     * Creates a folder with the given name and a log file that will be written to by the root logger.
     */
    private String initializeOutputFolderAndLogger(String anOutputFolderName) throws NullPointerException {
        String tmpOutputFolderPath = (new File(GlycosylationStatisticsTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator
                + anOutputFolderName + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + "Log.txt");
        } catch (IOException anIOException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.WARNING);
        return tmpOutputFolderPath;
    }

    /**
     * Creates and returns a PrintWriter instance that can be used to send output to a specified file.
     */
    private PrintWriter initializeOutputFile(String anOutputFolderPath, String anOutputFileName) throws IOException {
        String tmpOutputFilePath = anOutputFolderPath + anOutputFileName;
        File tmpOutputFile = new File(tmpOutputFilePath);
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        return new PrintWriter(tmpOutputFileWriter, true);
    }

    /**
     * Curation method for datasets of "synthetics". The structures in the given SMILES file are iterated and matched to
     * all COCONUT and ZINC "biogenic" subset molecules based on unique SMILES codes. All matching molecules are
     * filtered to obtain a dataset of "synthetic" molecules. Additionally, only one structure of every stereo-isomer group
     * in the original dataset is retained. The curated set is written to file as SMILES codes in the specified output folder.
     */
    private void filterZINCBiogenicAndCOCONUTMoleculesAndGroupStereoIsomersOfSMILESfile(
            String anOriginalDataSetFileName,
            int anOriginalDataSetSize,
            String aCuratedDataSetFileName,
            String anOutputFolderName)
            throws IllegalArgumentException, IOException {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        //Loading ZINC biogenic subset
        File tmpZincBiogenicSmilesFile = null;
        try {
            tmpZincBiogenicSmilesFile = new File(tmpClassLoader.getResource(GlycosylationStatisticsTest.ZINC_BIOGENIC_SUBSET_FILE_NAME).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("ZINC biogenic subset file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("ZINC biogenic subset file found at: " + tmpZincBiogenicSmilesFile.getAbsolutePath());
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(anOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        FileReader tmpZincBiogenicSmilesFileReader = new FileReader(tmpZincBiogenicSmilesFile);
        BufferedReader tmpZincBiogenicSmilesBufferedReader = new BufferedReader(tmpZincBiogenicSmilesFileReader);
        //key: unique SMILES string (does not encode stereochemistry!); object: ZINC id of the first stereo isomer encountered
        HashMap<String, String> tmpBiogenicSmilesMap = new HashMap<>((int)(308035*1.1), 1.0f); //308,035 molecules are in the dataset
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique); //Unique does not encode stereochemistry! Absolute would do
        String tmpZincBiogenicFileNextLine = "";
        String tmpZincBiogenicFileSmilesCode = "";
        String tmpID = "";
        int tmpZincBiogenicMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        System.out.println("Loading and processing ZINC biogenic subset now...");
        while (true) {
            try {
                tmpZincBiogenicFileNextLine = tmpZincBiogenicSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpZincBiogenicFileNextLine)) {
                    break;
                }
                if (tmpZincBiogenicFileNextLine.contains("SMILES")) {
                    continue;
                }
                tmpZincBiogenicMoleculesCounter++;
                if ((tmpZincBiogenicMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpZincBiogenicMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpZincBiogenicFileNextLine.split(" ");
                tmpZincBiogenicFileSmilesCode = tmpSmilesCodeAndId[0];
                tmpID = tmpSmilesCodeAndId[1];
                IAtomContainer tmpBiogenicMolecule = tmpSmiPar.parseSmiles(tmpZincBiogenicFileSmilesCode);
                String tmpCDKSmilesCode = tmpSmiGen.create(tmpBiogenicMolecule);
                if (!tmpBiogenicSmilesMap.containsKey(tmpCDKSmilesCode)) {
                    tmpBiogenicSmilesMap.put(tmpCDKSmilesCode, tmpID);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpZincBiogenicSmilesFileReader.close();
        tmpZincBiogenicSmilesBufferedReader.close();
        System.out.println("Processing of the biogenic subset done.");
        tmpOutputWriter.println("Processing of the biogenic subset done.");
        System.out.println(tmpZincBiogenicMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpZincBiogenicMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpBiogenicSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpBiogenicSmilesMap.size() + " SMILES codes have been put into memory.");
        System.out.println((tmpZincBiogenicMoleculesCounter - tmpBiogenicSmilesMap.size()) + " molecules were filtered because they are stereo-isomers of others.");
        tmpOutputWriter.println((tmpZincBiogenicMoleculesCounter - tmpBiogenicSmilesMap.size()) + " molecules were filtered because they are stereo-isomers of others.");
        //Loading COCONUT
        System.out.println("Loading and processing COCONUT now...");
        tmpOutputWriter.println("Loading and processing COCONUT now...");
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        //key: unique SMILES string (does not encode stereochemistry!); object: COCONUT id of the first stereo isomer encountered
        HashMap<String, String> tmpCoconutSmilesMap = new HashMap<>((int)(401624*1.1), 1.0f); //401,624 molecules are in the dataset
        Document tmpCurrentDoc;
        String tmpCoconutID = "";
        String tmpCoconutSmilesCode = "";
        int tmpCoconutMoleculesCounter = 0;
        tmpExceptionsCounter = 0;
        while (tmpCursor.hasNext()) {
            tmpCoconutID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                if (Objects.isNull(tmpCurrentDoc)) {
                    break;
                }
                tmpCoconutMoleculesCounter++;
                if ((tmpCoconutMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpCoconutMoleculesCounter + " molecules were processed already...");
                }
                tmpCoconutID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpCoconutSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                IAtomContainer tmpCoconutMolecule = tmpSmiPar.parseSmiles(tmpCoconutSmilesCode);
                String tmpCDKSmilesCode = tmpSmiGen.create(tmpCoconutMolecule);
                if (!tmpCoconutSmilesMap.containsKey(tmpCDKSmilesCode)) {
                    tmpCoconutSmilesMap.put(tmpCDKSmilesCode, tmpCoconutID);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpCoconutID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpCursor.close();
        System.out.println("Processing of COCONUT done.");
        tmpOutputWriter.println("Processing COCONUT done.");
        System.out.println(tmpCoconutMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpCoconutMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpCoconutSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpCoconutSmilesMap.size() + " SMILES codes have been put into memory.");
        System.out.println((tmpCoconutMoleculesCounter - tmpCoconutSmilesMap.size()) + " molecules were filtered because they are stereo-isomers of others.");
        tmpOutputWriter.println((tmpCoconutMoleculesCounter - tmpCoconutSmilesMap.size()) + " molecules were filtered because they are stereo-isomers of others.");
        //Loading analysed dataset
        System.out.println("Loading and processing the given dataset now...");
        tmpOutputWriter.println("Loading and processing the given dataset now...");
        File tmpRawSmilesFile = null;
        try {
            tmpRawSmilesFile = new File(tmpClassLoader.getResource(anOriginalDataSetFileName).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("Dataset file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Dataset found at: " + tmpRawSmilesFile.getAbsolutePath());
        FileReader tmpRawSmilesFileReader = new FileReader(tmpRawSmilesFile);
        BufferedReader tmpRawSmilesBufferedReader = new BufferedReader(tmpRawSmilesFileReader);
        //key: unique SMILES string (does not encode stereochemistry!)
        // object: HashMap of id of the first stereo isomer encountered and number of stereo-isomers encountered
        HashMap<String, HashMap<String, Object>> tmpSmilesMap = new HashMap<>((int)(anOriginalDataSetSize*1.1), 1.0f);
        String tmpRawFileNextLine = "";
        String tmpRawSmilesCode = "";
        tmpID = "";
        int tmpRawMoleculesCounter = 0;
        tmpExceptionsCounter = 0;
        int tmpMatchedBiogenicOrCOCONUTCounter = 0;
        while (true) {
            try {
                tmpRawFileNextLine = tmpRawSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpRawFileNextLine)) {
                    break;
                }
                if (tmpRawFileNextLine.contains("SMILES")) {
                    continue;
                }
                tmpRawMoleculesCounter++;
                if ((tmpRawMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpRawMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpRawFileNextLine.split(" ");
                tmpRawSmilesCode = tmpSmilesCodeAndId[0].trim();
                tmpID = tmpSmilesCodeAndId[1].trim();
                IAtomContainer tmpRawMoleculeFromSmiles = tmpSmiPar.parseSmiles(tmpRawSmilesCode);
                String tmpCDKSMILESCode = tmpSmiGen.create(tmpRawMoleculeFromSmiles);
                boolean tmpIsInBiogenic = tmpBiogenicSmilesMap.containsKey(tmpCDKSMILESCode);
                boolean tmpIsInCOCONUT = tmpCoconutSmilesMap.containsKey(tmpCDKSMILESCode);
                if (!(tmpIsInBiogenic || tmpIsInCOCONUT)) {
                    if (tmpSmilesMap.containsKey(tmpCDKSMILESCode)) {
                        HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpCDKSMILESCode);
                        int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
                        tmpInnerMap.put("FREQUENCY", (tmpStereoIsomerNumber + 1));
                        continue;
                    } else {
                        HashMap<String, Object> tmpInnerMap = new HashMap<>(3, 1.0f);
                        tmpInnerMap.put("ID", tmpID);
                        tmpInnerMap.put("FREQUENCY", 1);
                        tmpSmilesMap.put(tmpCDKSMILESCode, tmpInnerMap);
                        continue;
                    }
                } else {
                    tmpMatchedBiogenicOrCOCONUTCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpRawSmilesFileReader.close();
        tmpRawSmilesBufferedReader.close();
        System.out.println("Processing of the given dataset done.");
        tmpOutputWriter.println("Processing of the given dataset done.");
        System.out.println(tmpRawMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpRawMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        System.out.println(tmpMatchedBiogenicOrCOCONUTCounter + " molecules were filtered because they matched in biogenic or COCONUT.");
        tmpOutputWriter.println(tmpMatchedBiogenicOrCOCONUTCounter + " molecules were filtered because they matched in biogenic or COCONUT.");
        System.out.println((tmpRawMoleculesCounter - tmpMatchedBiogenicOrCOCONUTCounter - tmpSmilesMap.size())
                + " molecules were filtered because they were stereo isomers of already processed molecules.");
        tmpOutputWriter.println((tmpRawMoleculesCounter - tmpMatchedBiogenicOrCOCONUTCounter - tmpSmilesMap.size())
                + " molecules were filtered because they were stereo isomers of already processed molecules.");
        //Writing to file
        System.out.println("Writing the curated data set to file now...");
        PrintWriter tmpMoleculesWriter = this.initializeOutputFile(tmpOutputFolderPath, aCuratedDataSetFileName);
        int tmpCounter = 0;
        for (String tmpSmilesCode : tmpSmilesMap.keySet()) {
            if ((tmpCounter % 10000) == 0) {
                System.out.println(tmpCounter + " lines were processed already...");
            }
            HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpSmilesCode);
            int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
            tmpID = (String) tmpInnerMap.get("ID");
            tmpMoleculesWriter.println(tmpSmilesCode + " " + tmpID + " " + tmpStereoIsomerNumber);
            tmpCounter++;
        }
        System.out.println("Done, shutting down.");
        tmpMoleculesWriter.flush();
        tmpMoleculesWriter.close();
        tmpOutputWriter.flush();
        tmpOutputWriter.close();
    }

    /**
     * Curation method for datasets that should only be "flattened", i.e. stereo-isomers grouped. The structures in the
     * given SMILES file are iterated and only one structure of every stereo-isomer group
     * in the original dataset is written to the curated dataset. The curated set is written to file as SMILES codes in
     * the specified output folder.
     */
    private void groupStereoisomersOfSMILESfile(
            String anOriginalDataSetFileName,
            int anOriginalDataSetSize,
            String aCuratedDataSetFileName,
            String anOutputFolderName
    ) throws IllegalArgumentException, IOException {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(anOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //Loading analysed dataset
        System.out.println("Loading and processing the given dataset now...");
        tmpOutputWriter.println("Loading and processing the given dataset now...");
        File tmpOriginalDataSetSmilesFile = null;
        try {
            tmpOriginalDataSetSmilesFile = new File(tmpClassLoader.getResource(anOriginalDataSetFileName).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("Original data set file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Data set found at: " + tmpOriginalDataSetSmilesFile.getAbsolutePath());
        FileReader tmpOriginalDataSetSmilesFileReader = new FileReader(tmpOriginalDataSetSmilesFile);
        BufferedReader tmpOriginalDataSetSmilesBufferedReader = new BufferedReader(tmpOriginalDataSetSmilesFileReader);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique); //Unique does not encode stereochemistry! Absolute would do
        //key: unique SMILES string (does not encode stereochemistry!)
        //object: HashMap of ID of the first stereo isomer encountered and number of stereo-isomers encountered
        HashMap<String, HashMap<String, Object>> tmpSmilesMap = new HashMap<>((int)(anOriginalDataSetSize*1.1), 1.0f);
        String tmpOriginalDataSetFileNextLine = "";
        String tmpOriginalDataSetSmilesCode = "";
        String tmpID = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (true) {
            try {
                tmpOriginalDataSetFileNextLine = tmpOriginalDataSetSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpOriginalDataSetFileNextLine)) {
                    break;
                }
                if (tmpOriginalDataSetFileNextLine.contains("SMILES")) {
                    continue;
                }
                tmpMoleculesCounter++;
                if ((tmpMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpOriginalDataSetFileNextLine.split(" ");
                tmpOriginalDataSetSmilesCode = tmpSmilesCodeAndId[0].trim();
                tmpID = tmpSmilesCodeAndId[1].trim();
                IAtomContainer tmpMoleculeFromSmiles = tmpSmiPar.parseSmiles(tmpOriginalDataSetSmilesCode);
                String tmpCDKSMILESCode = tmpSmiGen.create(tmpMoleculeFromSmiles);
                if (tmpSmilesMap.containsKey(tmpCDKSMILESCode)) {
                    HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpCDKSMILESCode);
                    int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
                    tmpInnerMap.put("FREQUENCY", (tmpStereoIsomerNumber + 1));
                    continue;
                } else {
                    HashMap<String, Object> tmpInnerMap = new HashMap<>(3, 1.0f);
                    tmpInnerMap.put("ID", tmpID);
                    tmpInnerMap.put("FREQUENCY", 1);
                    tmpSmilesMap.put(tmpCDKSMILESCode, tmpInnerMap);
                    continue;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        tmpOriginalDataSetSmilesFileReader.close();
        tmpOriginalDataSetSmilesBufferedReader.close();
        System.out.println("Processing of the given dataset done.");
        tmpOutputWriter.println("Processing of the given dataset done.");
        System.out.println(tmpMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        tmpOutputWriter.println(tmpSmilesMap.size() + " SMILES codes have been put into memory.");
        //Writing to file
        System.out.println("Writing the curated data set to file now...");
        PrintWriter tmpMoleculesWriter = this.initializeOutputFile(tmpOutputFolderPath, aCuratedDataSetFileName);
        int tmpCounter = 0;
        tmpID = "";
        for (String tmpSmilesCode : tmpSmilesMap.keySet()) {
            if ((tmpCounter % 10000) == 0) {
                System.out.println(tmpCounter + " lines were processed already...");
            }
            HashMap<String, Object> tmpInnerMap = tmpSmilesMap.get(tmpSmilesCode);
            int tmpStereoIsomerNumber = (int) tmpInnerMap.get("FREQUENCY");
            tmpID = (String) tmpInnerMap.get("ID");
            tmpMoleculesWriter.println(tmpSmilesCode + " " + tmpID + " " + tmpStereoIsomerNumber);
            tmpCounter++;
        }
        System.out.println("Done, shutting down.");
        tmpMoleculesWriter.flush();
        tmpMoleculesWriter.close();
        tmpOutputWriter.flush();
        tmpOutputWriter.close();
    }

    /**
     * Calculates general glycosylation statistics on a given SMILES file. Calculated stats: How many molecules contain
     * (terminal/non-terminal)(circular/linear/both) sugars and how many molecules consist only of sugar units.
     * All statistics are printed to console and also compiled in an output file created in the given directory.
     * Test is ignored, if specified dataset cannot be found.
     */
    private void doBasicGlycoStatsOnSMILESfile(String aDatasetFileName, String anOutputFolderName) throws IllegalArgumentException, IOException {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSmilesFile = null;
        try {
            tmpSmilesFile = new File(tmpClassLoader.getResource(aDatasetFileName).getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("Data set could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println("Found input file at " + tmpSmilesFile.getAbsolutePath());
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(anOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        PrintWriter tmpSugarContainingMoleculesWriter = this.initializeOutputFile(tmpOutputFolderPath, "sugar-containing_molecules.txt");
        FileReader tmpSmilesFileReader = new FileReader(tmpSmilesFile);
        BufferedReader tmpSmilesBufferedReader = new BufferedReader(tmpSmilesFileReader);
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        String tmpFileNextLine = "";
        String tmpFileSmilesCode;
        String tmpID = "";
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsIDs = new ArrayList<>(10000);
        int tmpHasNoSugarsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsIDs = new ArrayList<>(7500);
        int tmpHasTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalCircularSugarsCNPs = new ArrayList<>(7500);
        int tmpHasNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasNonTerminalCircularSugarsCNPs = new ArrayList<>(7500);
        int tmpHasTerminalAndNonTerminalCircularSugarsCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsIDs = new ArrayList<>(2600);
        int tmpHasTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalLinearSugarsCNPs = new ArrayList<>(2600);
        int tmpHasNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasNonTerminalLinearSugarsCNPs = new ArrayList<>(2600);
        int tmpHasTerminalAndNonTerminalLinearSugarsCounter = 0;
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsIDs = new ArrayList<>(100);
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarIDs = new ArrayList<>(2000);
        int tmpTotalNrOfStereoIsomersSRUPos = 0;
        int tmpTotalNrOfStereoIsomersSRUNeg = 0;
        while (true) {
            try {
                tmpFileNextLine = tmpSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpFileNextLine)) {
                    break;
                }
                if (tmpFileNextLine.contains("SMILES")) {
                    continue;
                }
                if ((tmpMoleculesCounter % 10000) == 0) {
                    System.out.println(tmpMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpFileNextLine.split(" ");
                tmpFileSmilesCode = tmpSmilesCodeAndId[0];
                tmpID = tmpSmilesCodeAndId[1];
                tmpMoleculesCounter++;
                //use this value if the frequency is not given n the file
                int tmpNrOfStereoIsomers = 0;
                if (tmpSmilesCodeAndId.length == 3) {
                    tmpNrOfStereoIsomers = Integer.valueOf(tmpSmilesCodeAndId[2]);
                }
                tmpMolecule = tmpSmiPar.parseSmiles(tmpFileSmilesCode);
                if (!ConnectivityChecker.isConnected(tmpMolecule)) {
                    tmpMolecule = SugarRemovalUtility.selectBiggestUnconnectedFragment(tmpMolecule);
                }
                //note: per default, circular sugars having too few exocyclic oxygen atoms attached are not detected/removed/counted!
                // note also: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                boolean tmpHasAnyTypeOfSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyCircularSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY);
                boolean tmpHasAnyLinearSugar = tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY);
                if (tmpHasAnyTypeOfSugar) {
                    tmpHasAnyTypeOfSugarsCounter++;
                    tmpHasAnyTypeOfSugarsIDs.add(tmpID);
                    tmpTotalNrOfStereoIsomersSRUPos += tmpNrOfStereoIsomers;

                    if (tmpHasAnyCircularSugar) {
                        tmpHasCircularSugarsCounter++;
                        tmpHasCircularSugarsIDs.add(tmpID);
                        //terminal and non-terminal
                        List<IAtomContainer> tmpCircularSugarCandidatesList = tmpSugarRemovalUtil.getCircularSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfCircularSugarMoieties = tmpCircularSugarCandidatesList.size();
                        //note: circular moieties that become terminal after removal of a linear moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalCircularSugarMoieties = tmpSugarRemovalUtil.removeAndReturnCircularSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        int tmpNumberOfTerminalCircularSugarMoieties = tmpRemovedTerminalCircularSugarMoieties.size() - 1 ;
                        int tmpNumberOfNonTerminalCircularSugarMoieties = tmpNumberOfCircularSugarMoieties - tmpNumberOfTerminalCircularSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalCircularSugarMoieties >= 0);
                        if (tmpNumberOfTerminalCircularSugarMoieties > 0) {
                            tmpHasTerminalCircularSugarsCounter++;
                            tmpHasTerminalCircularSugarsCNPs.add(tmpID);
                        }
                        if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                            tmpHasNonTerminalCircularSugarsCounter++;
                            tmpHasNonTerminalCircularSugarsCNPs.add(tmpID);
                        }
                        if (tmpNumberOfTerminalCircularSugarMoieties > 0 && tmpNumberOfNonTerminalCircularSugarMoieties > 0) {
                            tmpHasTerminalAndNonTerminalCircularSugarsCounter++;
                        }
                    }
                    if (tmpHasAnyLinearSugar) {
                        tmpHasLinearSugarsCounter++;
                        tmpHasLinearSugarsIDs.add(tmpID);
                        //terminal and non-terminal
                        List<IAtomContainer> tmpLinearSugarCandidatesList = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpMolecule.clone());
                        int tmpNumberOfLinearSugarMoieties = tmpLinearSugarCandidatesList.size();
                        //note:linear moieties that become terminal after removal of a circular moiety are not counted here!
                        List<IAtomContainer> tmpRemovedTerminalLinearMoieties = tmpSugarRemovalUtil.removeAndReturnLinearSugars(tmpMolecule, true);
                        //-1 for the deglycosylated core at the beginning of the list
                        int tmpNumberOfTerminalLinearSugarMoieties = tmpRemovedTerminalLinearMoieties.size() - 1 ;
                        int tmpNumberOfNonTerminalLinearSugarMoieties = tmpNumberOfLinearSugarMoieties - tmpNumberOfTerminalLinearSugarMoieties;
                        Assert.assertTrue(tmpNumberOfNonTerminalLinearSugarMoieties >= 0);
                        if (tmpNumberOfTerminalLinearSugarMoieties > 0) {
                            tmpHasTerminalLinearSugarsCounter++;
                            tmpHasTerminalLinearSugarsCNPs.add(tmpID);
                        }
                        if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                            tmpHasNonTerminalLinearSugarsCounter++;
                            tmpHasNonTerminalLinearSugarsCNPs.add(tmpID);
                        }
                        if (tmpNumberOfTerminalLinearSugarMoieties > 0 && tmpNumberOfNonTerminalLinearSugarMoieties > 0) {
                            tmpHasTerminalAndNonTerminalLinearSugarsCounter++;
                        }
                    }
                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsIDs.add(tmpID);
                    }
                    if (tmpDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                        tmpBasicallyASugarIDs.add(tmpID);
                    }
                    String tmpSmilesCode = tmpSmiGen.create(tmpMolecule);
                    tmpSugarContainingMoleculesWriter.println(tmpSmilesCode + " " + tmpID);
                } else {
                    tmpHasNoSugarsCounter++;
                    tmpTotalNrOfStereoIsomersSRUNeg += tmpNrOfStereoIsomers;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        tmpOutputWriter.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        tmpOutputWriter.println("Sugar-containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        double tmpPercentage = ((double) tmpHasAnyTypeOfSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain sugars.");
        System.out.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        tmpOutputWriter.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpOutputWriter.println("Circular-sugar-containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        System.out.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        System.out.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        tmpOutputWriter.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        System.out.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpOutputWriter.println("Linear-sugar-containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        tmpOutputWriter.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        System.out.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        System.out.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        tmpOutputWriter.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        tmpOutputWriter.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpOutputWriter.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        System.out.println("Total number of stereo isomers of molecules without sugars in the original dataset: "
                + tmpTotalNrOfStereoIsomersSRUNeg);
        tmpOutputWriter.println("Total number of stereo isomers of molecules without sugars in the original dataset: "
                + tmpTotalNrOfStereoIsomersSRUNeg);
        System.out.println("Stereo isomers per molecule on average: "
                + ((double) tmpTotalNrOfStereoIsomersSRUNeg / (double) tmpHasNoSugarsCounter));
        tmpOutputWriter.println("Stereo isomers per molecule on average: "
                + ((double) tmpTotalNrOfStereoIsomersSRUNeg / (double) tmpHasNoSugarsCounter));
        System.out.println("Total number of stereo isomers of molecules with sugars in the original dataset: "
                + tmpTotalNrOfStereoIsomersSRUPos);
        tmpOutputWriter.println("Total number of stereo isomers of molecules with sugars in the original dataset: "
                + tmpTotalNrOfStereoIsomersSRUPos);
        System.out.println("Stereo isomers per molecule on average: "
                + ((double) tmpTotalNrOfStereoIsomersSRUPos / (double) tmpHasAnyTypeOfSugarsCounter));
        tmpOutputWriter.println("Stereo isomers per molecule on average: "
                + ((double) tmpTotalNrOfStereoIsomersSRUPos / (double) tmpHasAnyTypeOfSugarsCounter));
        tmpOutputWriter.println();
        tmpOutputWriter.println("Sugar-containing molecules: " + tmpHasAnyTypeOfSugarsIDs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Circular-sugar-containing molecules: " + tmpHasCircularSugarsIDs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal circular sugar containing molecules: " + tmpHasTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Non-terminal circular sugar containing molecules: " + tmpHasNonTerminalCircularSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Linear-sugar-containing molecules: " + tmpHasLinearSugarsIDs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Terminal linear sugars containing molecules: " + tmpHasTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Non-terminal linear sugar containing molecules: " + tmpHasNonTerminalLinearSugarsCNPs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules containing both circular and linear sugars: " + tmpHasCircularAndLinearSugarsIDs);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Basically a sugar: " + tmpBasicallyASugarIDs);
        tmpOutputWriter.flush();
        tmpSugarContainingMoleculesWriter.flush();
        tmpSmilesFileReader.close();
        tmpSmilesBufferedReader.close();
        tmpOutputWriter.close();
        tmpSugarContainingMoleculesWriter.close();
    }
    //</editor-fold>
}

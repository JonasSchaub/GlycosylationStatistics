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

/**
 * TODO:
 * - write doc
 * - update ZINC
 * - test whether the sugar-containing molecules in ZINC are NPs or are actually also part of COCONUT (really, how?)
 * - Add more stats, any ideas??
 * - subdivide the detected linear sugars in rings somehow, this number is odd! By size? By size of the rings they are part of?
 * - include only needed CDK modules
 */

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
//TODO
import com.mongodb.client.*;
import de.unijena.cheminf.deglycosylation.SugarRemovalUtility;
import org.bson.Document;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.hash.HashGeneratorMaker;
import org.openscience.cdk.hash.MoleculeHashGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.ringsearch.RingSearch;
import org.openscience.cdk.smiles.SmiFlavor;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

//TODO
import java.awt.*;
import java.io.*;
import java.math.RoundingMode;
import java.text.NumberFormat;
//TODO
import java.util.*;
import java.util.List;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;

/**
 * TODO
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
     * Name of the MongoDB database
     */
    private static final String DATABASE_NAME = "COCONUT2020september09";

    /**
     * Collection from the database to load
     */
    private static final String COLLECTION_NAME = "uniqueNaturalProduct";

    /**
     * TODO
     */
    private static final String SDF_NAME = "COCONUT_DB_september_09.sdf";

    /**
     * Name of the output folder
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
     * TODO
     */
    private static final String OUTPUT_FILE_SEPARATOR = ";";

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(GlycosylationStatisticsTest.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Tests involving databases">
    //<editor-fold desc="COCONUT">
    /**
     * TODO
     * Ignored because many files are written
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
                    tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_deglycosylated.png");
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
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
     * TODO
     * Consistency tested and approved
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
        //List<String> tmpHasNoSugarsCNPs = new ArrayList<>(400000);
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList<>(3000);
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList<>(2000);
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
                    //tmpHasNoSugarsCNPs.add(tmpID);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
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
        //tmpOutputWriter.println("No sugar containing molecules: " + tmpHasNoSugarsCNPs);
        //tmpOutputWriter.println();
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
     * TODO
     * Consistency tested and approved
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
        List<String> tmpHasTerminalCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasNonTerminalCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasOnlyTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasOnlyNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasTerminalAndNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalCircularSugarsCNPs = new ArrayList<>(5000);
        int tmpHasOnlyCircularSugarsCounter = 0;
        List<String> tmpHasOnlyCircularSugarsCNPs = new ArrayList<>(50000);
        int tmpHasGlycosidicBondCounter = 0;
        List<String> tmpHasGlycosidicBondCNPs = new ArrayList<>(50000);
        int tmpHasGlycosidicBondOnTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnTerminalSugarCNPs = new ArrayList<>(50000);
        int tmpHasGlycosidicBondOnNonTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnNonTerminalSugarCNPs = new ArrayList<>(50000);
        int tmpGlycosidicBondExemptionCounter = 0;
        List<String> tmpGlycosidicBondExemptionCNPs = new ArrayList<>(120);
        int tmpHasSpiroSugarCounter = 0;
        List<String> tmpHasSpiroSugarCNPs = new ArrayList<>(100);
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
                //continue;
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
        tmpOutputWriter.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasOnlyCircularSugarsCounter) + "molecules " +
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i));
            tmpCSVWriter.println(i + ":" + tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugars += tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i);
        } */
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpCSVheavyAtomCountWriter.println(i + ":" + tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars1 += tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i);
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpCSVcarbonAtomCountWriter.println(i + ":" + tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars2 += tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i);
        } */
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
     * TODO: Depict molecules with linear sugars in cycles highlighted.
     */
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
        //TODO
        PrintWriter tmpCSVWriter = this.initializeOutputFile(tmpOutputFolderPath, "LinSugarsCarbonAtomCountFrequencies.csv");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //tmpSugarRemovalUtil.setDetectSpiroRingsAsCircularSugarsSetting(true);
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
        int tmpLinearSugarMoietiesInRingsCounter2 = 0;
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
                    //so far not needed here
                    //tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    int tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidatesExcludingCycles.size();
                    int tmpNumberOfLinearSugarsInCycles = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    Assert.assertTrue(tmpNumberOfLinearSugarsInCycles >= 0);
                    //this will be the list of candidates that only get detected if cyclic atoms are included
                    if (tmpNumberOfLinearSugarsInCycles > 0) {
                        //System.out.println(tmpNumberOfLinearSugarsInCycles);
                        tmpLinearSugarMoietiesInRingsCounter += tmpNumberOfLinearSugarsInCycles;
                        List<IAtomContainer> tmpLinearCandidatesActuallyInRings = new ArrayList<>(tmpNumberOfLinearSugarsInCycles * 2);
                        int[][] tmpAdjList = GraphUtil.toAdjList(tmpMolecule);
                        RingSearch tmpRingSearch = new RingSearch(tmpMolecule, tmpAdjList);
                        for (IAtomContainer tmpCandidate : tmpLinearCandidatesIncludingCycles) {
                            for (IAtom tmpAtom : tmpCandidate.atoms()) {
                                if (tmpRingSearch.cyclic(tmpAtom)) {
                                    tmpLinearCandidatesActuallyInRings.add(tmpCandidate);
                                    //move on with the next candidate
                                    break;
                                }
                            }
                        }
                        /*tmpDepictionGenerator.withHighlight(tmpLinearCandidatesActuallyInRings, Color.BLUE)
                                .withSize(2000, 2000)
                                .withFillToFit()
                                .depict(tmpMolecule)
                                .writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");*/
                        //System.out.println(tmpLinearCandidatesActuallyInRings.size());
                        /*for (int i = 0; i < tmpLinearCandidatesActuallyInRings.size(); i++) {
                            IAtomContainer tmpCandidate = tmpLinearCandidatesActuallyInRings.get(i);
                            int tmpCarbonAtomCountTotal = 0;
                            int tmpCarbonAtomCountInRing = 0;
                            for (IAtom tmpAtom : tmpCandidate.atoms()) {
                                boolean tmpIsCarbon = tmpAtom.getSymbol().equals("C");
                                boolean tmpIsCyclic = tmpRingSearch.cyclic(tmpAtom);
                                if (tmpIsCarbon) {
                                    tmpCarbonAtomCountTotal++;
                                    if (tmpIsCyclic) {
                                        tmpCarbonAtomCountInRing++;
                                    }
                                }
                            }
                            int tmpExocyclicCarbonAtomCount = tmpCarbonAtomCountTotal - tmpCarbonAtomCountInRing;
                            if (tmpExocyclicCarbonAtomCount < tmpSugarRemovalUtil.getLinearSugarCandidateMinSizeSetting()) {
                                tmpLinearCandidatesActuallyInRings.remove(i);
                                i = i - 1;
                            }
                        }
                        //System.out.println(tmpLinearCandidatesActuallyInRings.size());
                        //System.out.println();*/
                        tmpLinearSugarMoietiesInRingsCounter2 += tmpLinearCandidatesActuallyInRings.size();
                        //Assert.assertTrue(tmpLinearCandidatesActuallyInRings.size() == tmpNumberOfLinearSugarsInCycles);
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
                //continue;
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
        System.out.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter2);
        tmpOutputWriter.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter2);
        System.out.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinSugInRingsLostInRemovalOfCircSugCounter);
        tmpOutputWriter.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinSugInRingsLostInRemovalOfCircSugCounter);
        tmpOutputWriter.println();
        tmpOutputWriter.println("Molecules that lost a linear sugar in a ring after removal of circular sugars: "
                + tmpLinSugInRingsLostInRemovalOfCircSugCNPs);
        tmpOutputWriter.flush();
        tmpCSVWriter.flush();
        tmpCursor.close();
        tmpOutputWriter.close();
        tmpCSVWriter.close();
    }

    /**
     * TODO
     * Consistency tested and approved
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
                //continue;
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpHowManyMoleculesHaveHowManySugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i));
            tmpCSVmoietyNrFreqWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i));
            tmpTotalOfSugarContainingMolecules += tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i);
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i));
            tmpCSVcircMoietyNrFreqWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i);
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.keySet()) {
            System.out.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i));
            tmpCSVcircMoietyGlyBondNrFreqWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i));
            tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i);
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ":" + tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i));
            tmpCSVlinMoietyNrFreqWriter.println(i + ":" + tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i);
        } */
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
        /* //Not occurring ratios are missing in the printout:
        List<String> tmpSortedKeySet = new ArrayList<>(tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.keySet());
        Collections.sort(tmpSortedKeySet);
        for (String i : tmpSortedKeySet) {
            System.out.println(i + ":" + tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(i));
            tmpCSVExoCycOxRatioFreqWriter.println(i + ":" + tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(i));
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i));
            tmpCSVExoCycOxFuranosesFreqWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i));
        }*/
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i));
            tmpCSVExoCycOxPyranosesFreqWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i));
        } */
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
        /* //Not occurring ratios are missing in the printout:
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.keySet()) {
            System.out.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i));
            tmpOutputWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i));
            tmpCSVExoCycOxHeptosesFreqWriter.println(i + ":" + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i));
        } */
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
     * TODO
     * Consistency tested and approved
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
                //continue;
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
     * TODO
     *
     * @throws Exception if anything goes wrong
     */
    //@Ignore
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
                //continue;
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
        System.out.println("Linear sugar patter SMILES codes and detected frequencies: ");
        tmpOutputWriter.println("Linear sugar patter SMILES codes and detected frequencies: ");
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
     * TODO
     * Note: Numbers give appearance of moiety, not how many molecules have this moiety!
     * Note: Per default, circular sugars having too few exocyclic oxygen atoms attached are not counted!
     * Problem: These moieties are only the rings without any further information!
     * @throws Exception
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
        String tmpCSVmoietyFreqFileHeader = "hash" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "frequency" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "firstOrigin";
        tmpCSVmoietyFreqWriter.println(tmpCSVmoietyFreqFileHeader);
        tmpCSVmoietyFreqWriter.flush();
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        //maybe adjust this
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        MoleculeHashGenerator tmpHashGenerator = new HashGeneratorMaker().depth(16)
                .elemental()
                .charged()
                .molecular();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpOutput = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasCircularSugarsCounter = 0;
        int tmpDifferentMoietiesCounter = 0;
        HashMap<Long, HashMap> tmpCircularSugarMoietiesMap = new HashMap<>(2000, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                //using the clean smiles without hydrogens here because of the depiction and the hashing of molecules
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");//tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
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
                        long tmpHashCode = tmpHashGenerator.generate(tmpMoiety);
                        String tmpMoietySmilesCode;
                        try {
                            tmpMoietySmilesCode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException anException) {
                            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                            tmpMoietySmilesCode = "[exception]";
                        }
                        if (tmpCircularSugarMoietiesMap.containsKey(tmpHashCode)) {
                            HashMap<String, Object> tmpInnerMap = tmpCircularSugarMoietiesMap.get(tmpHashCode);
                            int tmpFrequency = (int)tmpInnerMap.get("FREQUENCY");
                            tmpInnerMap.put("FREQUENCY", tmpFrequency + 1);
                        } else {
                            tmpDifferentMoietiesCounter++;
                            HashMap<String, Object> tmpInnerMap = new HashMap<>(5,1);
                            tmpInnerMap.put("FREQUENCY", 1);
                            tmpInnerMap.put("FIRST_ORIGIN", tmpID);
                            tmpInnerMap.put("SMILES", tmpMoietySmilesCode);
                            tmpCircularSugarMoietiesMap.put(tmpHashCode, tmpInnerMap);
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
                //continue;
            }
        }
        Long[] tmpKeySet = tmpCircularSugarMoietiesMap.keySet().toArray(new Long[0]);
        Arrays.sort(tmpKeySet, new Comparator<Long>() {
            public int compare(Long aFirstKey, Long aSecondKey) {
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
        for (long tmpHashCode : tmpKeySet) {
            String tmpEntry = tmpHashCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR;
            HashMap<String, Object> tmpInnerMap = tmpCircularSugarMoietiesMap.get(tmpHashCode);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("SMILES") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("FREQUENCY") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat((String)tmpInnerMap.get("FIRST_ORIGIN"));
            tmpCSVmoietyFreqWriter.println(tmpEntry);
            if ((int)tmpInnerMap.get("FREQUENCY") > 9) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpInnerMap.get("SMILES")))
                        .writeTo(tmpOutputFolderPath + File.separator + tmpInnerMap.get("FREQUENCY")
                                + "_" + tmpHashCode +".png");
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
     * TODO
     * Note: Numbers give appearance of moiety, not how many molecules have this moiety!
     * Note: per default, linear sugars in rings, those too small or too big, and acidic linear sugars are not detected/removed/counted!
     * To consider: Additonal functionalities on the sugars are not reflected in these moieties.
     * @throws Exception
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
        String tmpCSVmoietyFreqFileHeader = "hash" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "SMILES" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "frequency" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR
                + "firstOrigin";
        tmpCSVmoietyFreqWriter.println(tmpCSVmoietyFreqFileHeader);
        tmpCSVmoietyFreqWriter.flush();
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        //maybe adjust this
        tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SmilesGenerator tmpSmiGen = new SmilesGenerator(SmiFlavor.Unique);
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
        MoleculeHashGenerator tmpHashGenerator = new HashGeneratorMaker().depth(16)
                .elemental()
                .charged()
                .molecular();
        Document tmpCurrentDoc;
        String tmpID;
        String tmpSmilesCode;
        IAtomContainer tmpMolecule;
        IAtomContainer tmpMoleculeWithoutSugars;
        String tmpOutput = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasLinearSugarsCounter = 0;
        int tmpDifferentMoietiesCounter = 0;
        HashMap<Long, HashMap> tmpLinearSugarMoietiesMap = new HashMap<>(2000, 0.9f);
        while (tmpCursor.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculesCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                //using the clean smiles without hydrogens here because of the depiction and the hashing of molecules
                tmpSmilesCode = tmpCurrentDoc.getString("clean_smiles");//tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
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
                        long tmpHashCode = tmpHashGenerator.generate(tmpMoiety);
                        String tmpMoietySmilesCode;
                        try {
                            tmpMoietySmilesCode = tmpSmiGen.create(tmpMoiety);
                        } catch (CDKException anException) {
                            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                            tmpMoietySmilesCode = "[exception]";
                        }
                        if (tmpLinearSugarMoietiesMap.containsKey(tmpHashCode)) {
                            HashMap<String, Object> tmpInnerMap = tmpLinearSugarMoietiesMap.get(tmpHashCode);
                            int tmpFrequency = (int)tmpInnerMap.get("FREQUENCY");
                            tmpInnerMap.put("FREQUENCY", tmpFrequency + 1);
                        } else {
                            tmpDifferentMoietiesCounter++;
                            HashMap<String, Object> tmpInnerMap = new HashMap<>(5,1);
                            tmpInnerMap.put("FREQUENCY", 1);
                            tmpInnerMap.put("FIRST_ORIGIN", tmpID);
                            tmpInnerMap.put("SMILES", tmpMoietySmilesCode);
                            tmpLinearSugarMoietiesMap.put(tmpHashCode, tmpInnerMap);
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
                //continue;
            }
        }
        Long[] tmpKeySet = tmpLinearSugarMoietiesMap.keySet().toArray(new Long[0]);
        Arrays.sort(tmpKeySet, new Comparator<Long>() {
            public int compare(Long aFirstKey, Long aSecondKey) {
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
        for (long tmpHashCode : tmpKeySet) {
            String tmpEntry = tmpHashCode + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR;
            HashMap<String, Object> tmpInnerMap = tmpLinearSugarMoietiesMap.get(tmpHashCode);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("SMILES") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat(tmpInnerMap.get("FREQUENCY") + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR);
            tmpEntry = tmpEntry.concat((String)tmpInnerMap.get("FIRST_ORIGIN"));
            tmpCSVmoietyFreqWriter.println(tmpEntry);
            if ((int)tmpInnerMap.get("FREQUENCY") > 9) {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpInnerMap.get("SMILES")))
                        .writeTo(tmpOutputFolderPath + File.separator + tmpInnerMap.get("FREQUENCY")
                                + "_" + tmpHashCode +".png");
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
     * TODO
     * Note: Ignoring stereo-chemistry and grouping sugar moieties that are the same without stereo-chemistry
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
        tmpCSVSRUPositiveWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
        tmpCSVSRUNegativeWriter.println("ID" + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + "Frequency");
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
        MoleculeHashGenerator tmpHashGenerator = new HashGeneratorMaker().depth(16)
                .elemental()
                .charged()
                .molecular();
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), true);
        HashMap<Long, HashMap<String, Object>> tmpSRUPositiveSugarPatterns = new HashMap<>(344, 1);
        HashMap<Long, HashMap<String, Object>> tmpSRUNegativeSugarPatterns = new HashMap<>(344, 1);
        String tmpReviewSugarID;
        int tmpReviewSugarsCounter = 0;
        int tmpReviewDataExceptionsCounter = 0;
        while (tmpReader.hasNext()) {
            tmpReviewSugarID = "[unidentified]";
            try {
                IAtomContainer tmpReviewSugar = tmpReader.next();
                tmpReviewSugarsCounter++;
                tmpReviewSugarID = tmpReviewSugar.getProperty("Name");
                long tmpHashCode = tmpHashGenerator.generate(tmpReviewSugar);
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
                }
                tmpMap.put("SMILES", tmpSmilesCode);
                //TODO: Also test whether it gets completely removed?
                boolean tmpHasSugars = tmpSugarRemovalUtil.hasCircularSugars(tmpReviewSugar);
                if (tmpHasSugars) {
                    if (tmpSRUPositiveSugarPatterns.containsKey(tmpHashCode)) {
                        HashMap tmpInnerMap = tmpSRUPositiveSugarPatterns.get(tmpHashCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUPositiveSugarPatterns.put(tmpHashCode, tmpMap);
                    }
                } else {
                    if (tmpSRUNegativeSugarPatterns.containsKey(tmpHashCode)) {
                        HashMap tmpInnerMap = tmpSRUNegativeSugarPatterns.get(tmpHashCode);
                        tmpInnerMap.put("ID", tmpInnerMap.get("ID") + "_" + tmpReviewSugarID);
                    } else {
                        tmpSRUNegativeSugarPatterns.put(tmpHashCode, tmpMap);
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpReviewSugarID, anException);
                tmpReviewDataExceptionsCounter++;
                //continue;
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
                //continue;
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
        HashMap<String, Object>[] tmpSRUPositivePatternsMapArray = tmpSRUPositiveSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUPositivePatternsMapArray, new Comparator<HashMap>() {
            public int compare(HashMap aFirstMap, HashMap aSecondMap) {
                int tmpFirstFrequency = (int)aFirstMap.get("FREQUENCY");
                int tmpSecondFrequency = (int)aSecondMap.get("FREQUENCY");
                return Integer.compare(tmpFirstFrequency, tmpSecondFrequency) * (-1);
            }
        });
        HashMap<String, Object>[] tmpSRUNegativePatternsMapArray = tmpSRUNegativeSugarPatterns.values().toArray(new HashMap[0]);
        Arrays.sort(tmpSRUNegativePatternsMapArray, new Comparator<HashMap>() {
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
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUPositiveWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && (String)tmpReviewSugarMap.get("SMILES") != "[generation_failed]") {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpReviewSugarMap.get("SMILES")))
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
            int tmpFrequency = (int) tmpReviewSugarMap.get("FREQUENCY");
            tmpCSVSRUNegativeWriter.println(tmpName + GlycosylationStatisticsTest.OUTPUT_FILE_SEPARATOR + tmpFrequency);
            if (tmpFrequency > 9 && (String)tmpReviewSugarMap.get("SMILES") != "[generation_failed]") {
                tmpDepictionGenerator.withSize(2000, 2000)
                        .withFillToFit()
                        .depict(tmpSmiPar.parseSmiles((String)tmpReviewSugarMap.get("SMILES")))
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
     * TODO
     *
     * @throws Exception if anything goes wrong
     */
    @Ignore
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
        //List<String> tmpHasNoSugarsCNPs = new ArrayList<>(400000);
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
                    //tmpHasNoSugarsCNPs.add(tmpID);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
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
        //tmpOutputWriter.println("No sugar containing molecules: " + tmpHasNoSugarsCNPs);
        //tmpOutputWriter.println();
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
    //</editor-fold>

    //<editor-fold desc="ZINC">
    /**
     * TODO
     *
     * @throws Exception if anything goes wrong
     */
    @Ignore
    @Test
    public void zincStatsBasicsTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpZincSDFile = null;
        try {
            tmpZincSDFile = new File(tmpClassLoader.getResource("zinc-all-for-sale.sdf").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("ZINC SDF could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpZincSDFile.getAbsolutePath());
        final String tmpSpecificOutputFolderName = "zinc_stats_basics_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpZincSDFile), DefaultChemObjectBuilder.getInstance(), true);
        String tmpID;
        IAtomContainer tmpMolecule;
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsCNPs = new ArrayList<>(10000);
        int tmpHasNoSugarsCounter = 0;
        //List<String> tmpHasNoSugarsCNPs = new ArrayList<>(3100000);
        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList<>(7500);
        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList<>(2600);
        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList<>(100);
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarCNPs = new ArrayList<>(2000);
        while (tmpReader.hasNext()) {
            tmpID = "[unidentified]";
            try {
                tmpMolecule = tmpReader.next();
                tmpMoleculesCounter++;
                tmpID = tmpMolecule.getProperty("zinc_id");
                tmpMolecule.setTitle(tmpID);
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
                    //tmpHasNoSugarsCNPs.add(tmpID);
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
            }
        }
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        tmpOutputWriter.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
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
        //tmpOutputWriter.println("No sugar containing molecules: " + tmpHasNoSugarsCNPs);
        //tmpOutputWriter.println();
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
     * TODO
     */
    @Test
    public void zincDatasetCurationTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpZincBiogenicSmilesFile = null;
        try {
            tmpZincBiogenicSmilesFile = new File(tmpClassLoader.getResource("ZINC_biogenic_subset_2020_Okt_16.txt").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("ZINC biogenic subset file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpZincBiogenicSmilesFile.getAbsolutePath());
        final String tmpSpecificOutputFolderName = "zinc_curation_test";
        //Prints output folder to console
        String tmpOutputFolderPath = this.initializeOutputFolderAndLogger(tmpSpecificOutputFolderName);
        PrintWriter tmpOutputWriter = this.initializeOutputFile(tmpOutputFolderPath, "Output.txt");
        FileReader tmpZincBiogenicSmilesFileReader = new FileReader(tmpZincBiogenicSmilesFile);
        BufferedReader tmpZincBiogenicSmilesBufferedReader = new BufferedReader(tmpZincBiogenicSmilesFileReader);
        HashSet<String> tmpBiogenicSmilesSet = new HashSet<>(135400, 1); //135,335 molecules are in the dataset
        String tmpNextLine = "";
        String tmpSmilesCode;
        String tmpID = "";
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        while (true) {
            try {
                tmpNextLine = tmpZincBiogenicSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpNextLine)) {
                    break;
                }
                tmpMoleculesCounter++;
                if ((tmpMoleculesCounter % 1000) == 0) {
                    System.out.println(tmpMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpNextLine.split(" ");
                tmpSmilesCode = tmpSmilesCodeAndId[0];
                tmpID = tmpSmilesCodeAndId[1];
                tmpBiogenicSmilesSet.add(tmpSmilesCode);
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
            }
        }
        tmpZincBiogenicSmilesFileReader.close();
        tmpZincBiogenicSmilesBufferedReader.close();
        System.out.println("Processing of the biogenic subset done.");
        System.out.println(tmpMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println("Loading and processing for-sale subset now...");
        tmpOutputWriter.println("Processing of the biogenic subset done.");
        tmpOutputWriter.println(tmpMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println("Loading and processing for-sale subset now...");
        File tmpZincForSaleSmilesFile = null;
        try {
            tmpZincForSaleSmilesFile = new File(tmpClassLoader.getResource("ZINC_for-sale_subset_2020_Okt_16.txt").getFile());
        } catch (NullPointerException aNullPointerException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aNullPointerException.toString(), aNullPointerException);
            System.out.println("ZINC for sale subset file could not be found. Test is ignored.");
            Assume.assumeTrue(false);
        }
        System.out.println(tmpZincForSaleSmilesFile.getAbsolutePath());
        FileReader tmpZincForSaleSmilesFileReader = new FileReader(tmpZincForSaleSmilesFile);
        BufferedReader tmpZincForSaleSmilesBufferedReader = new BufferedReader(tmpZincForSaleSmilesFileReader);
        File tmpSdFile = new File(tmpOutputFolderPath + File.separator + "ZINC_for-sale_without_biogenics.sdf");
        FileWriter tmpSdfFileWriter = new FileWriter(tmpSdFile);
        BufferedWriter tmpSDFBufferedWriter = new BufferedWriter(tmpSdfFileWriter);
        SDFWriter tmpSdfWriter = new SDFWriter(tmpSDFBufferedWriter);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        tmpNextLine = "";
        tmpSmilesCode = "";
        tmpID = "";
        tmpMoleculesCounter = 0;
        tmpExceptionsCounter = 0;
        int tmpMoleculesWrittenToSDFCounter = 0;
        while (true) {
            try {
                tmpNextLine = tmpZincForSaleSmilesBufferedReader.readLine();
                if (Objects.isNull(tmpNextLine)) {
                    break;
                }
                tmpMoleculesCounter++;
                if ((tmpMoleculesCounter % 1000) == 0) {
                    System.out.println(tmpMoleculesCounter + " lines were processed already...");
                }
                String[] tmpSmilesCodeAndId = tmpNextLine.split(" ");
                tmpSmilesCode = tmpSmilesCodeAndId[0];
                tmpID = tmpSmilesCodeAndId[1];
                boolean tmpIsBiogenic = tmpBiogenicSmilesSet.contains(tmpSmilesCode);
                if (!tmpIsBiogenic) {
                    IAtomContainer tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                    tmpSdfWriter.write(tmpMolecule);
                    tmpMoleculesWrittenToSDFCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                //continue;
            }
        }
        tmpZincForSaleSmilesFileReader.close();
        tmpZincForSaleSmilesBufferedReader.close();
        tmpSdfWriter.close();
        tmpSDFBufferedWriter.flush();
        tmpSDFBufferedWriter.close();
        tmpSdfFileWriter.flush();
        tmpSdfFileWriter.close();
        System.out.println("Processing of the for-sale subset done.");
        System.out.println(tmpMoleculesCounter + " molecules were processed.");
        System.out.println(tmpExceptionsCounter + " exceptions occurred.");
        System.out.println(tmpMoleculesWrittenToSDFCounter + " molecules were written to the SD file.");
        tmpOutputWriter.println("Processing of the for-sale subset done.");
        tmpOutputWriter.println(tmpMoleculesCounter + " molecules were processed.");
        tmpOutputWriter.println(tmpExceptionsCounter + " exceptions occurred.");
        tmpOutputWriter.println(tmpMoleculesWrittenToSDFCounter + " molecules were written to the SD file.");
    }
    //</editor-fold>
    //</editor-fold>
    //
    //<editor-fold desc="Private methods">
    /**
     * TODO
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
     * TODO
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
     * TODO
     */
    private PrintWriter initializeOutputFile(String anOutputFolderPath, String anOutputFileName) throws IOException {
        String tmpOutputFilePath = anOutputFolderPath + anOutputFileName;
        File tmpOutputFile = new File(tmpOutputFilePath);
        FileWriter tmpOutputFileWriter = new FileWriter(tmpOutputFile);
        return new PrintWriter(tmpOutputFileWriter, true);
    }
    //</editor-fold>
}

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
 * - Split up Coconut stats test
 * - Redo statistics with new COCONUT version and new algo (also ZINC)
 * - test whether the sugar-containing molecules in ZINC are NPs or are actually also part of COCONUT
 * - Do explicit imports
 * - write doc
 */

import com.mongodb.MongoClientSettings;
import com.mongodb.MongoTimeoutException;
import com.mongodb.ServerAddress;
import com.mongodb.client.*;
import de.unijena.cheminf.deglycosylation.SugarRemovalUtility;
import org.bson.Document;
import org.junit.Assert;
import org.junit.Assume;
import org.junit.Ignore;
import org.junit.Test;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.isomorphism.DfPattern;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.NumberFormat;
import java.util.*;
import java.util.logging.FileHandler;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.logging.SimpleFormatter;
import java.util.stream.Collectors;

/**
 * TODO
 *
 * @author Jonas Schaub
 * @version 1.0.0.0
 */
public class GlycosylationStatisticsTest extends SugarRemovalUtility {
    //<editor-fold defaultstate="collapsed" desc="Private static final constants">
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
    private static final String DATABASE_NAME = "COCONUTaugust17";

    /**
     * Collection from the database to load
     */
    private static final String COLLECTION_NAME = "uniqueNaturalProduct";

    /**
     * Name of the output folder
     */
    private static final String OUTPUT_FOLDER_NAME = "GlycosylationStatisticsTest_Output";

    /**
     * Name of the document variable that contains an ID of the given molecule, "coconut_id" or "inchikey"
     */
    private static final String ID_KEY = "coconut_id";

    /**
     * Name of the document variable that contains the SMILES code string of the given molecule, "clean_smiles" or "smiles"
     */
    private static final String SMILES_CODE_KEY = "smiles";

    /**
     * Logger of this class
     */
    private static final Logger LOGGER = Logger.getLogger(GlycosylationStatisticsTest.class.getName());
    //</editor-fold>
    //
    //<editor-fold desc="Tests involving databases">
    /**
     * TODO
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void visualInspectionOfCoconutTest() throws Exception {
        final String tmpSpecificOutputFolderName = "coconut_inspection_test";
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        String tmpOutputFolderPath = (new File(GlycosylationStatisticsTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator
                + tmpSpecificOutputFolderName + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        DepictionGenerator tmpDepictionGenerator = new DepictionGenerator();
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
        //All settings in default
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
                tmpID = tmpCurrentDoc.getString(GlycosylationStatisticsTest.ID_KEY);
                tmpSmilesCode = tmpCurrentDoc.getString(GlycosylationStatisticsTest.SMILES_CODE_KEY);
                tmpMolecule = tmpSmiPar.parseSmiles(tmpSmilesCode);
                tmpMolecule.setTitle(tmpID);
                IAtomContainer tmpDeglycosylatedClone = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    tmpDepictionGenerator.depict(tmpMolecule).writeTo(tmpOutputFolderPath + File.separator + tmpID + ".png");
                    tmpDepictionGenerator.depict(tmpDeglycosylatedClone).writeTo(tmpOutputFolderPath + File.separator + tmpID + "_deglycosylated.png");
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpDeglycosylatedClone.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsLinearSugarsCounter++;
                    }
                    if (tmpDeglycosylatedClone.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                    }
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpSugarContainingMoleculesCounter);
        System.out.println("Linear-sugar-containing molecules counter: " + tmpContainsLinearSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpContainsCircularSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpCursor.close();
    }

    /**
     * TODO split this up into multiple tests!
     */
    @Ignore
    @Test
    public void CoconutStatsTest() throws Exception {

        //<editor-fold desc="Connection to MongoDB">
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        //</editor-fold>

        //<editor-fold desc="Output folder and logging setup">
        final String tmpSpecificOutputFolderName = "coconut_stats_test";
        String tmpOutputFolderPath = (new File(GlycosylationStatisticsTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator
                + tmpSpecificOutputFolderName + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
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
        //</editor-fold>

        NumberFormat tmpRatioOutputFormat = NumberFormat.getInstance(Locale.US);
        tmpRatioOutputFormat.setMaximumFractionDigits(1);
        tmpRatioOutputFormat.setRoundingMode(RoundingMode.DOWN);

        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);

        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;

        //<editor-fold desc="Counter and map definitions">
        int tmpMoleculesCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpHasAnyTypeOfSugarsCounter = 0;
        List<String> tmpHasAnyTypeOfSugarsCNPs = new ArrayList<>(50000);
        int tmpHasNoSugarsCounter = 0;
        //List<String> tmpHasNoSugarsCNPs = new ArrayList<>(400000);

        int tmpHasCircularSugarsCounter = 0;
        List<String> tmpHasCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasNonTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasOnlyTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasOnlyNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasTerminalAndNonTerminalCircularSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalCircularSugarsCNPs = new ArrayList(5000);
        int tmpHasOnlyCircularSugarsCounter = 0;
        List<String> tmpHasOnlyCircularSugarsCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondCounter = 0;
        List<String> tmpHasGlycosidicBondCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondOnTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnTerminalSugarCNPs = new ArrayList(50000);
        int tmpHasGlycosidicBondOnNonTerminalSugarCounter = 0;
        List<String> tmpHasGlycosidicBondOnNonTerminalSugarCNPs = new ArrayList(50000);
        int tmpGlycosidicBondExemptionCounter = 0;
        List<String> tmpGlycosidicBondExemptionCNPs = new ArrayList(120);

        int tmpHasLinearSugarsCounter = 0;
        List<String> tmpHasLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasNonTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasOnlyTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasOnlyNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasOnlyNonTerminalLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasTerminalAndNonTerminalLinearSugarsCounter = 0;
        List<String> tmpHasTerminalAndNonTerminalLinearSugarsCNPs = new ArrayList(10);
        int tmpHasOnlyLinearSugarsCounter = 0;
        List<String> tmpHasOnlyLinearSugarsCNPs = new ArrayList(3000);
        int tmpHasLinearSugarsInRingCounter = 0;
        List<String> tmpHasLinearSugarsInRingCNPs = new ArrayList(500);

        int tmpHasCircularAndLinearSugarsCounter = 0;
        List<String> tmpHasCircularAndLinearSugarsCNPs = new ArrayList(2000);

        //number of molecules that are basically sugars, circular or linear, polymer or single unit
        int tmpBasicallyASugarCounter = 0;
        List<String> tmpBasicallyASugarCNPs = new ArrayList(2000);
        int tmpBasicallyACircularSugarCounter = 0;
        List<String> tmpBasicallyACircularSugarCNPs = new ArrayList(2000);
        int tmpBasicallyALinearSugarCounter = 0;
        List<String> tmpBasicallyALinearSugarCNPs = new ArrayList(2000);
        int tmpBasicallyASingleSugarUnitCounter = 0; //circular or linear
        List<String> tmpBasicallyASingleSugarUnitCNPs = new ArrayList(1000);
        int tmpBasicallyASingleCircularSugarCounter = 0;
        List<String> tmpBasicallyASingleCircularSugarCNPs = new ArrayList(1000);
        int tmpBasicallyASingleLinearSugarCounter = 0;
        List<String> tmpBasicallyASingleLinearSugarCNPs = new ArrayList(1000);

        int tmpCircularSugarMoietiesCounter = 0;
        int tmpTerminalCircularSugarMoietiesCounter = 0;
        int tmpNonTerminalCircularSugarMoietiesCounter = 0;
        int tmpCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest have no glycosidic bond
        int tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter = 0; //the rest are non-terminal

        int tmpLinearSugarMoietiesCounter = 0;
        int tmpTerminalLinearSugarMoietiesCounter = 0;
        int tmpNonTerminalLinearSugarMoietiesCounter = 0;
        int tmpLinearSugarMoietiesInRingsCounter = 0;

        HashMap<Integer, Integer> tmpFrequenciesOfSizesOfCircularSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap = new HashMap(10, 0.9f);

        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManySugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap = new HashMap(10, 0.9f);

        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<Integer, Integer> tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap = new HashMap(10, 0.9f);
        HashMap<String, Integer> tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap = new HashMap(10, 0.9f);
        int tmpUnexpectedRingSizeCounter = 0;

        int tmpLinearSugarsDetectedInCircularSugarsCounter = 0;
        //</editor-fold>

        //<editor-fold desc="Iteration of molecules">
        while (tmpCursor.hasNext()) {
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
                    tmpHasAnyTypeOfSugarsCNPs.add(tmpID);

                    int tmpNumberOfCircularAndLinearSugarMoieties = tmpSugarRemovalUtil.getNumberOfCircularAndLinearSugars(tmpMolecule);
                    if (!tmpHowManyMoleculesHaveHowManySugarMoietiesMap.containsKey(tmpNumberOfCircularAndLinearSugarMoieties)) {
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, 1);
                    } else {
                        Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(tmpNumberOfCircularAndLinearSugarMoieties);
                        tmpHowManyMoleculesHaveHowManySugarMoietiesMap.put(tmpNumberOfCircularAndLinearSugarMoieties, tmpCurrentListValue + 1);
                    }

                    if (tmpHasAnyCircularSugar) {

                        //<editor-fold desc="Analysis of molecules having circular sugars">
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

                        if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.containsKey(tmpNumberOfCircularSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(tmpNumberOfCircularSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.put(tmpNumberOfCircularSugarMoieties, tmpCurrentValue + 1);
                        }

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
                            tmpTerminalCircularSugarMoietiesCounter += tmpNumberOfTerminalCircularSugarMoieties;

                            if (tmpNumberOfNonTerminalCircularSugarMoieties == 0) {
                                tmpHasOnlyTerminalCircularSugarsCounter++;
                                tmpHasOnlyTerminalCircularSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfNonTerminalCircularSugarMoieties > 0) {

                            tmpHasNonTerminalCircularSugarsCounter++;
                            tmpHasNonTerminalCircularSugarsCNPs.add(tmpID);
                            tmpNonTerminalCircularSugarMoietiesCounter += tmpNumberOfNonTerminalCircularSugarMoieties;

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

                            //the first is the overall counter, the second one is specific for this molecule
                            tmpCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfGlycosidicBonds;

                            if (!tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.containsKey(tmpNumberOfGlycosidicBonds)) {
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, 1);
                            } else {
                                Integer tmpCurrentListValue = tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(tmpNumberOfGlycosidicBonds);
                                tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.put(tmpNumberOfGlycosidicBonds, tmpCurrentListValue + 1);
                            }

                            if (tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond > 0) {
                                tmpHasGlycosidicBondOnTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnTerminalSugarCNPs.add(tmpID);
                                //the first is the overall counter, the second one is specific for this molecule
                                tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter += tmpNumberOfTerminalSugarMoietiesWithGlycosidicBond;
                            }

                            if (tmpNumberOfNonTerminalSugarMoietiesWithGlycosidicBond > 0) {
                                tmpHasGlycosidicBondOnNonTerminalSugarCounter++;
                                tmpHasGlycosidicBondOnNonTerminalSugarCNPs.add(tmpID);
                            }

                        }
                        //</editor-fold>

                    }

                    //<editor-fold desc="Analysis of molecules having linear sugars">
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

                        //the first is the overall counter, the second one is specific for this molecule
                        tmpLinearSugarMoietiesCounter += tmpNumberOfLinearSugarMoieties;

                        if (!tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.containsKey(tmpNumberOfLinearSugarMoieties)) {
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, 1);
                        } else {
                            Integer tmpCurrentValue = tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(tmpNumberOfLinearSugarMoieties);
                            tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.put(tmpNumberOfLinearSugarMoieties, tmpCurrentValue + 1);
                        }

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

                            tmpHasTerminalLinearSugarsCounter++;
                            tmpHasTerminalLinearSugarsCNPs.add(tmpID);
                            tmpTerminalLinearSugarMoietiesCounter += tmpNumberOfTerminalLinearSugarMoieties;

                            if (tmpNumberOfNonTerminalLinearSugarMoieties == 0) {
                                tmpHasOnlyTerminalLinearSugarsCounter++;
                                tmpHasOnlyTerminalLinearSugarsCNPs.add(tmpID);
                            }

                        }

                        if (tmpNumberOfNonTerminalLinearSugarMoieties > 0) {

                            tmpHasNonTerminalLinearSugarsCounter++;
                            tmpHasNonTerminalLinearSugarsCNPs.add(tmpID);
                            tmpNonTerminalLinearSugarMoietiesCounter += tmpNumberOfNonTerminalLinearSugarMoieties;

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
                    //</editor-fold>

                    if (tmpHasAnyCircularSugar && tmpHasAnyLinearSugar) {
                        tmpHasCircularAndLinearSugarsCounter++;
                        tmpHasCircularAndLinearSugarsCNPs.add(tmpID);
                    }

                    //<editor-fold desc="Analysis of molecules that are basically sugars">
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
                    //</editor-fold>

                } else {
                    tmpHasNoSugarsCounter++;
                    //tmpHasNoSugarsCNPs.add(tmpID);
                }

                //<editor-fold desc="Analysis of all possible sugar rings and their exocyclic oxygen counts">
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
                //</editor-fold>

                //<editor-fold desc="Analysis of linear sugars in rings">
                //leaving default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                IAtomContainer tmpNewClone = tmpMolecule.clone();
                List<IAtomContainer> tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                int tmpListSizeWithCandidatesInCycles = tmpLinearCandidates.size();
                if (tmpListSizeWithCandidatesInCycles > 0) {
                    //this.removeSugarCandidatesWithCyclicAtoms(tmpLinearCandidates, tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    int tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidates.size();
                    int tmpNumberOfLinearSugarsInCycles = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    if (tmpNumberOfLinearSugarsInCycles > 0) {
                        tmpHasLinearSugarsInRingCounter++;
                        tmpHasLinearSugarsInRingCNPs.add(tmpID);
                        tmpLinearSugarMoietiesInRingsCounter += tmpNumberOfLinearSugarsInCycles;
                    }
                    //leaving default further!
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(false);
                    tmpSugarRemovalUtil.removeCircularSugars(tmpNewClone, false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpListSizeWithCandidatesInCycles = tmpLinearCandidates.size();
                    //this.removeSugarCandidatesWithCyclicAtoms(tmpLinearCandidates, tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                    tmpLinearCandidates = tmpSugarRemovalUtil.getLinearSugarCandidates(tmpNewClone);
                    tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(true);
                    tmpListSizeWithoutCandidatesInCycles = tmpLinearCandidates.size();
                    int tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars = tmpListSizeWithCandidatesInCycles - tmpListSizeWithoutCandidatesInCycles;
                    int tmpLinearSugarsThatWereDetectedInCircularSugars = tmpNumberOfLinearSugarsInCycles - tmpNumberOfLinearSugarsInCyclesWithoutCircularSugars;
                    tmpLinearSugarsDetectedInCircularSugarsCounter += tmpLinearSugarsThatWereDetectedInCircularSugars;
                    //back to this default
                    tmpSugarRemovalUtil.setRemoveOnlyTerminalSugarsSetting(true);
                }
                //back to default settings!
                tmpSugarRemovalUtil.setDetectLinearSugarsInRingsSetting(false);
                //</editor-fold>

            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
            }

        } //end of while()
        //</editor-fold>

        //<editor-fold desc="Printout">
        System.out.println();
        System.out.println("Done.");
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println();
        System.out.println("Molecules counter: " + tmpMoleculesCounter);
        System.out.println();
        System.out.println("Sugar containing molecules counter: " + tmpHasAnyTypeOfSugarsCounter);
        double tmpPercentage = ((double) tmpHasAnyTypeOfSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain sugars.");
        System.out.println("No sugar containing molecules counter: " + tmpHasNoSugarsCounter);
        System.out.println();
        System.out.println("Circular sugar containing molecules counter: " + tmpHasCircularSugarsCounter);
        tmpPercentage = ((double) tmpHasCircularSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain circular sugars.");
        System.out.println("Only circular sugar containing molecules counter: " + tmpHasOnlyCircularSugarsCounter);
        System.out.println("The rest contain linear sugars only or both types of sugars (see below)");
        System.out.println("Terminal circular sugars containing molecules counter: " + tmpHasTerminalCircularSugarsCounter);
        System.out.println("Only terminal circular sugar containing molecules counter: " + tmpHasOnlyTerminalCircularSugarsCounter);
        System.out.println("Non-terminal circular sugar containing molecules counter: " + tmpHasNonTerminalCircularSugarsCounter);
        System.out.println("Only non-terminal circular sugar containing molecules counter: " + tmpHasOnlyNonTerminalCircularSugarsCounter);
        System.out.println("Terminal and non-terminal circular sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalCircularSugarsCounter);
        System.out.println("Circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondCounter);
        System.out.println("The remaining " + (tmpHasCircularSugarsCounter - tmpHasGlycosidicBondCounter) + " molecules " +
                "only have circular sugar moieties that are NOT attached via a glycosidic bond.");
        System.out.println("Terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on non-terminal circular sugar moieties.");
        System.out.println("Non-terminal circular sugar with glycosidic bond containing molecules counter: " + tmpHasGlycosidicBondOnNonTerminalSugarCounter);
        System.out.println("The remaining " + (tmpHasGlycosidicBondCounter - tmpHasGlycosidicBondOnNonTerminalSugarCounter) + " molecules " +
                "only have glycosidic bonds on terminal circular sugar moieties.");
        int tmpHasBoth = tmpHasGlycosidicBondOnNonTerminalSugarCounter + tmpHasGlycosidicBondOnTerminalSugarCounter - tmpHasGlycosidicBondCounter;
        System.out.println(tmpHasBoth + " molecules have both, terminal and non-terminal sugar moieties attached via a glycosidic bond.");
        System.out.println("Molecules that qualify for the glycosidic bond exemption counter: " + tmpGlycosidicBondExemptionCounter);
        System.out.println();
        System.out.println("Detected circular sugar moieties counter: " + tmpCircularSugarMoietiesCounter);
        System.out.println("Detected terminal circular sugar moieties counter: " + tmpTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected non-terminal circular sugar moieties counter: " + tmpNonTerminalCircularSugarMoietiesCounter);
        System.out.println("Detected circular sugar moieties that have a glycosidic bond counter: " + tmpCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesCounter - tmpCircularSugarMoietiesWithGlycosidicBondCounter) + " circular sugar moieties do not have a glycosidic bond.");
        System.out.println("Detected circular sugar moieties that have a glycosidic bond and are terminal counter: " + tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter);
        System.out.println((tmpCircularSugarMoietiesWithGlycosidicBondCounter - tmpTerminalCircularSugarMoietiesWithGlycosidicBondCounter)
                + " circular sugar moieties that have a glycosidic bond are non-terminal.");
        System.out.println();
        System.out.println("Linear sugar containing molecules counter: " + tmpHasLinearSugarsCounter);
        tmpPercentage = ((double) tmpHasLinearSugarsCounter / (double) tmpMoleculesCounter) * 100;
        System.out.println(tmpPercentage + " % of molecules contain linear sugars.");
        System.out.println("Only linear sugar containing molecules counter: " + tmpHasOnlyLinearSugarsCounter);
        System.out.println("Terminal linear sugars containing molecules counter: " + tmpHasTerminalLinearSugarsCounter);
        System.out.println("Only terminal linear sugar containing molecules counter: " + tmpHasOnlyTerminalLinearSugarsCounter);
        System.out.println("Non-terminal linear sugar containing molecules counter: " + tmpHasNonTerminalLinearSugarsCounter);
        System.out.println("Only non-terminal linear sugar containing molecules counter: " + tmpHasOnlyNonTerminalLinearSugarsCounter);
        System.out.println("Terminal and non-terminal linear sugar containing molecules counter: " + tmpHasTerminalAndNonTerminalLinearSugarsCounter);
        System.out.println("Linear sugar moieties in rings containing molecules counter: " + tmpHasLinearSugarsInRingCounter);
        System.out.println();
        System.out.println("Detected linear sugar moieties counter: " + tmpLinearSugarMoietiesCounter);
        System.out.println("Detected terminal linear sugar moieties counter: " + tmpTerminalLinearSugarMoietiesCounter);
        System.out.println("Detected non-terminal linear sugar moieties counter: " + tmpNonTerminalLinearSugarMoietiesCounter);
        System.out.println("Detected linear sugar moieties that are part of rings counter: " + tmpLinearSugarMoietiesInRingsCounter);
        System.out.println();
        System.out.println("Molecules containing both circular and linear sugars counter: " + tmpHasCircularAndLinearSugarsCounter);
        System.out.println();
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        System.out.println("Basically a single sugar unit counter: " + tmpBasicallyASingleSugarUnitCounter);
        System.out.println("Basically a circular sugar counter: " + tmpBasicallyACircularSugarCounter);
        System.out.println("Basically a single circular sugar counter: " + tmpBasicallyASingleCircularSugarCounter);
        System.out.println("Basically a linear sugar counter: " + tmpBasicallyALinearSugarCounter);
        System.out.println("Basically a single linear sugar counter: " + tmpBasicallyASingleLinearSugarCounter);
        System.out.println();
        System.out.println("How many molecules have how many sugars: ");
        int tmpTotalOfSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManySugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i));
            tmpTotalOfSugarContainingMolecules += tmpHowManyMoleculesHaveHowManySugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many circular sugars: ");
        int tmpTotalOfCircularSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many circular sugars attached via a glycosidic bond: ");
        int tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i));
            tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules += tmpHowManyMoleculesHaveHowManyCircularSugarMoietiesWithGlycosidicBondMap.get(i);
        }
        System.out.println();
        System.out.println("How many molecules have how many linear sugars: ");
        int tmpTotalOfLinearSugarContainingMolecules = 0;
        for (int i : tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugarContainingMolecules += tmpHowManyMoleculesHaveHowManyLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= heavy atom count) frequency distribution of circular sugars: ");
        int tmpTotalOfCircularSugars = 0;
        for (int i : tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i));
            tmpTotalOfCircularSugars += tmpFrequenciesOfSizesOfCircularSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= heavy atom count) frequency distribution of linear sugars: ");
        int tmpTotalOfLinearSugars1 = 0;
        for (int i : tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars1 += tmpFrequenciesOfHeavyAtomCountsOfLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Size (= carbon atom count) frequency distribution of linear sugars: ");
        int tmpTotalOfLinearSugars2 = 0;
        for (int i : tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i));
            tmpTotalOfLinearSugars2 += tmpFrequenciesOfCarbonAtomCountsOfLinearSugarMoietiesMap.get(i);
        }
        System.out.println();
        System.out.println("Terminal Sugars counter: " + (tmpTerminalCircularSugarMoietiesCounter + tmpTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpTerminalLinearSugarMoietiesCounter + " of these are linear");
        System.out.println();
        System.out.println("Non-terminal Sugars counter: " + (tmpNonTerminalCircularSugarMoietiesCounter + tmpNonTerminalLinearSugarMoietiesCounter));
        System.out.println(tmpNonTerminalCircularSugarMoietiesCounter + " of these are circular");
        System.out.println(tmpNonTerminalLinearSugarMoietiesCounter + " of these are linear");
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atoms to atoms in ring ratios of circular sugars: ");
        List<String> tmpSortedKeySet = tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.keySet().stream().collect(Collectors.toList());
        Collections.sort(tmpSortedKeySet);
        for (String i : tmpSortedKeySet) {
            System.out.println(i + ": " + tmpFrequenciesOfAttachedExocyclicOxygenAtomsRatiosMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 5-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn5MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 6-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn6MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Frequency distribution of exocyclic oxygen atom counts of 7-membered circular sugars: ");
        for (int i : tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.keySet()) {
            System.out.println(i + ": " + tmpFrequenciesOfNumbersOfAttachedExocyclicOxygenAtomsOn7MemberedRingsMap.get(i));
        }
        System.out.println();
        System.out.println("Number of circular sugar moieties that had an unexpected ring size (should be zero!): " + tmpUnexpectedRingSizeCounter);
        System.out.println();
        System.out.println("Number of detected linear sugars in rings that got lost through the removal of circular " +
                "sugars counter: " + tmpLinearSugarsDetectedInCircularSugarsCounter);
        System.out.println();

        //Example print out of CNP list elements:
        //System.out.println("These molecules are basically sugars:");
        //System.out.println(Arrays.toString(tmpBasicallyASingleSugarUnitCNPs.toArray()));
        //System.out.println(tmpBasicallyASingleSugarUnitCNPs); //shorter alternative
        //</editor-fold>

        tmpCursor.close();

        //<editor-fold desc="Tests for consistency">
        Assert.assertEquals(tmpMoleculesCounter, tmpHasNoSugarsCounter + tmpHasAnyTypeOfSugarsCounter);
        Assert.assertEquals(tmpHasAnyTypeOfSugarsCounter, tmpHasOnlyCircularSugarsCounter
                + tmpHasOnlyLinearSugarsCounter
                + tmpHasCircularAndLinearSugarsCounter);
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpHasTerminalAndNonTerminalCircularSugarsCounter
                + tmpHasOnlyTerminalCircularSugarsCounter
                + tmpHasOnlyNonTerminalCircularSugarsCounter);
        Assert.assertTrue(tmpHasBoth >= 0);
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpHasTerminalAndNonTerminalLinearSugarsCounter
                + tmpHasOnlyTerminalLinearSugarsCounter
                + tmpHasOnlyNonTerminalLinearSugarsCounter);
        Assert.assertEquals(tmpHasAnyTypeOfSugarsCounter, tmpTotalOfSugarContainingMolecules);
        Assert.assertEquals(tmpHasCircularSugarsCounter, tmpTotalOfCircularSugarContainingMolecules);
        Assert.assertEquals(tmpHasGlycosidicBondCounter, tmpTotalOfCircularSugarWithGlycosidicBondContainingMolecules);
        Assert.assertEquals(tmpHasLinearSugarsCounter, tmpTotalOfLinearSugarContainingMolecules);
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTerminalCircularSugarMoietiesCounter
                + tmpNonTerminalCircularSugarMoietiesCounter);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTerminalLinearSugarMoietiesCounter
                + tmpNonTerminalLinearSugarMoietiesCounter);
        Assert.assertEquals(tmpCircularSugarMoietiesCounter, tmpTotalOfCircularSugars);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars1);
        Assert.assertEquals(tmpLinearSugarMoietiesCounter, tmpTotalOfLinearSugars2);
        Assert.assertEquals(tmpBasicallyASingleSugarUnitCounter,
                tmpBasicallyASingleCircularSugarCounter + tmpBasicallyASingleLinearSugarCounter);
        //</editor-fold>
    }

    /**
     * TODO
     *
     * @throws Exception
     */
    @Ignore
    @Test
    public void linearSugarPatternsAppearanceTest() throws Exception {
        MongoCursor<Document> tmpCursor = null;
        try {
            //prints to console if connection was successful
            tmpCursor = this.getCOCONUTMongoCursorForIteration();
        } catch (MongoTimeoutException aMongoTimeoutException) {
            GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, aMongoTimeoutException.toString(), aMongoTimeoutException);
            System.out.println("Timed out while trying to connect to MongoDB. Test is ignored.");
            Assume.assumeTrue(false);
        }
        SmilesParser tmpSmiPar = new SmilesParser(DefaultChemObjectBuilder.getInstance());
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        //to also test their appearance
        tmpSugarRemovalUtil.setDetectLinearAcidicSugarsSetting(true);
        //Note: Here, additional molecules could be added to the list to also test them
        List<String> tmpLinearSugarsList = tmpSugarRemovalUtil.getLinearSugarPatternsList();
        List<List<Object>> tmpLinearSugarPatterns = new ArrayList<>(tmpLinearSugarsList.size());
        for (String tmpLinearSugarString : tmpLinearSugarsList) {
            List<Object> tmpList = new ArrayList<>(4);
            tmpList.add(0, tmpLinearSugarString);
            tmpList.add(1, DfPattern.findSubstructure(tmpSmiPar.parseSmiles(tmpLinearSugarString)));
            tmpList.add(2, 0);
            tmpLinearSugarPatterns.add(tmpList);
        }
        Document tmpCurrentDoc = null;
        String tmpID = "";
        String tmpSmilesCode = "";
        IAtomContainer tmpMolecule = null;
        int tmpMoleculeCounter = 0;
        int tmpContainsLinearSugarCounter = 0;
        while (tmpCursor.hasNext()) {
            try {
                tmpCurrentDoc = tmpCursor.next();
                tmpMoleculeCounter++;
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
                    tmpContainsLinearSugarCounter++;
                }
            } catch (Exception anException) {
                GlycosylationStatisticsTest.LOGGER.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Molecules that have at least one match counter: " + tmpContainsLinearSugarCounter);
        for (List<Object> tmpEntry : tmpLinearSugarPatterns) {
            System.out.println((String)tmpEntry.get(0) + " " + (int)tmpEntry.get(2));
        }
        tmpCursor.close();
    }

    /**
     * TODO
     */
    @Ignore
    @Test
    public void zincInspectionTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpZincSDFile = new File(tmpClassLoader.getResource("zinc-all-for-sale.sdf").getFile());
        System.out.println(tmpZincSDFile.getAbsolutePath());
        Logger tmpLogger = Logger.getLogger(GlycosylationStatisticsTest.class.getName());
        final String tmpSpecificOutputFolderName = "zinc_inspection_test";
        String tmpOutputFolderPath = (new File(GlycosylationStatisticsTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator
                + tmpSpecificOutputFolderName + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + "Log.txt");
        } catch (IOException anIOException) {
            tmpLogger.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.WARNING);
        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpZincSDFile), DefaultChemObjectBuilder.getInstance(), true);
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        IAtomContainer tmpMolecule = null;
        String tmpID = null;
        while (tmpReader.hasNext()) {
            try {
                tmpMolecule = tmpReader.next();
                tmpID = tmpMolecule.getProperty("zinc_id");
                tmpMoleculeCounter++;
                tmpMolecule = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsLinearSugarsCounter++;
                    }
                    if (tmpMolecule.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                    }
                }
            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpSugarContainingMoleculesCounter);
        System.out.println("Linear-sugar-containing molecules counter: " + tmpContainsLinearSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpContainsCircularSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpReader.close();
    }

    /**
     * TODO
     */
    @Ignore
    @Test
    public void coconutSdfTest() throws Exception {
        ClassLoader tmpClassLoader = this.getClass().getClassLoader();
        File tmpSDFile = new File(tmpClassLoader.getResource("COCONUT_DB_august_17.sdf").getFile());
        System.out.println(tmpSDFile.getAbsolutePath());
        Logger tmpLogger = Logger.getLogger(GlycosylationStatisticsTest.class.getName());
        final String tmpSpecificOutputFolderName = "coconut_sdf_test";
        String tmpOutputFolderPath = (new File(GlycosylationStatisticsTest.OUTPUT_FOLDER_NAME)).getAbsolutePath() + File.separator
                + tmpSpecificOutputFolderName + File.separator;
        File tmpOutputFolderFile = new File(tmpOutputFolderPath);
        if (!tmpOutputFolderFile.exists()) {
            tmpOutputFolderFile.mkdirs();
        }
        System.out.println("Output directory: " + tmpOutputFolderPath);
        FileHandler tmpLogFileHandler = null;
        try {
            tmpLogFileHandler = new FileHandler(tmpOutputFolderPath + "Log.txt");
        } catch (IOException anIOException) {
            tmpLogger.log(Level.SEVERE, anIOException.toString(), anIOException);
            System.out.println("An exception occurred while setting up the log file. Logging will be done in default configuration.");
        }
        tmpLogFileHandler.setLevel(Level.ALL);
        tmpLogFileHandler.setFormatter(new SimpleFormatter());
        Logger.getLogger("").addHandler(tmpLogFileHandler);
        Logger.getLogger("").setLevel(Level.WARNING);
        //Done for reproducibility
        SugarRemovalUtility tmpSugarRemovalUtil = new SugarRemovalUtility();
        tmpSugarRemovalUtil.setAddPropertyToSugarContainingMoleculesSetting(true);
        IteratingSDFReader tmpReader = new IteratingSDFReader(new FileInputStream(tmpSDFile), DefaultChemObjectBuilder.getInstance(), true);
        int tmpMoleculeCounter = 0;
        int tmpExceptionsCounter = 0;
        int tmpSugarContainingMoleculesCounter = 0;
        int tmpContainsLinearSugarsCounter = 0;
        int tmpContainsCircularSugarsCounter = 0;
        int tmpBasicallyASugarCounter = 0;
        IAtomContainer tmpMolecule = null;
        String tmpID = null;
        while (tmpReader.hasNext()) {
            try {
                tmpMolecule = tmpReader.next();
                tmpID = tmpMolecule.getProperty(GlycosylationStatisticsTest.ID_KEY);
                tmpMoleculeCounter++;
                tmpMolecule = tmpSugarRemovalUtil.removeCircularAndLinearSugars(tmpMolecule, true);
                if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_SUGAR_PROPERTY_KEY) == true) {
                    tmpSugarContainingMoleculesCounter++;
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_CIRCULAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsCircularSugarsCounter++;
                    }
                    if ((boolean)tmpMolecule.getProperty(SugarRemovalUtility.CONTAINS_LINEAR_SUGAR_PROPERTY_KEY) == true) {
                        tmpContainsLinearSugarsCounter++;
                    }
                    if (tmpMolecule.isEmpty()) {
                        tmpBasicallyASugarCounter++;
                    }
                }
            } catch (Exception anException) {
                tmpLogger.log(Level.SEVERE, anException.toString() + " ID: " + tmpID, anException);
                tmpExceptionsCounter++;
                continue;
            }
        }
        System.out.println("Done.");
        System.out.println("Molecules counter: " + tmpMoleculeCounter);
        System.out.println("Exceptions counter: " + tmpExceptionsCounter);
        System.out.println("Sugar-containing molecules counter: " + tmpSugarContainingMoleculesCounter);
        System.out.println("Linear-sugar-containing molecules counter: " + tmpContainsLinearSugarsCounter);
        System.out.println("Circular-sugar-containing molecules counter: " + tmpContainsCircularSugarsCounter);
        System.out.println("Basically a sugar counter: " + tmpBasicallyASugarCounter);
        tmpReader.close();
    }
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
        MongoCursor<Document> tmpCursor = null;
        tmpCursor = tmpCollection.find().iterator();
        System.out.println("Connection to MongoDB successful.");
        System.out.println("Collection " + tmpCollection.getNamespace().getCollectionName() + " in database " + tmpDatabase.getName() + " is loaded.");
        return tmpCursor;
    }
    //</editor-fold>
}

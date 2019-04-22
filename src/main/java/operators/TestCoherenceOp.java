package operators;

import com.bc.ceres.core.ProgressMonitor;
import org.esa.s1tbx.commons.Sentinel1Utils;
import org.esa.s1tbx.insar.gpf.InterferogramOp;
import org.esa.snap.core.datamodel.*;
import org.esa.snap.core.dataop.dem.ElevationModel;
import org.esa.snap.core.dataop.resamp.ResamplingFactory;
import org.esa.snap.core.gpf.Operator;
import org.esa.snap.core.gpf.OperatorException;
import org.esa.snap.core.gpf.Tile;
import org.esa.snap.core.gpf.annotations.OperatorMetadata;
import org.esa.snap.core.gpf.annotations.Parameter;
import org.esa.snap.core.gpf.annotations.SourceProduct;
import org.esa.snap.core.gpf.annotations.TargetProduct;
import org.esa.snap.core.util.ProductUtils;
import org.esa.snap.dem.dataio.DEMFactory;
import org.esa.snap.dem.dataio.FileElevationModel;
import org.esa.snap.engine_utilities.datamodel.AbstractMetadata;
import org.esa.snap.engine_utilities.datamodel.PosVector;
import org.esa.snap.engine_utilities.datamodel.Unit;
import org.esa.snap.engine_utilities.eo.Constants;
import org.esa.snap.engine_utilities.eo.GeoUtils;
import org.esa.snap.engine_utilities.gpf.InputProductValidator;
import org.esa.snap.engine_utilities.gpf.OperatorUtils;
import org.esa.snap.engine_utilities.gpf.StackUtils;
import org.esa.snap.engine_utilities.gpf.TileIndex;
import org.jblas.ComplexDouble;
import org.jblas.ComplexDoubleMatrix;
import org.jblas.DoubleMatrix;
import org.jblas.MatrixFunctions;
import org.jlinda.core.GeoPoint;
import org.jlinda.core.Orbit;
import org.jlinda.core.Point;
import org.jlinda.core.SLCImage;
import org.jlinda.core.geom.DemTile;
import org.jlinda.core.geom.TopoPhase;
import org.jlinda.core.utils.*;

import javax.media.jai.BorderExtender;
import java.awt.*;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;

@OperatorMetadata(alias = "Coherence",
        category = "Radar/Interferometric/Products",
        authors = "Petar Marinkovic, Jun Lu",
        version = "1.0",
        copyright = "Copyright (C) 2013 by PPO.labs",
        description = "Estimate coherence from stack of coregistered images")
public class TestCoherenceOp extends Operator {

    @SourceProduct
    public Product sourceProduct;

    @TargetProduct
    private Product targetProduct;

    @Parameter(interval = "(1, 90]",
            description = "Size of coherence estimation window in Azimuth direction",
            defaultValue = "10",
            label = "Coherence Azimuth Window Size")
    private int cohWinAz = 10;

    @Parameter(interval = "(1, 90]",
            description = "Size of coherence estimation window in Range direction",
            defaultValue = "10",
            label = "Coherence Range Window Size")
    private int cohWinRg = 10;

    @Parameter(defaultValue = "false", label = "Subtract flat-earth phase in coherence phase")
    private boolean subtractFlatEarthPhase = false;

    @Parameter(valueSet = {"1", "2", "3", "4", "5", "6", "7", "8"},
            description = "Order of 'Flat earth phase' polynomial",
            defaultValue = "5",
            label = "Degree of \"Flat Earth\" polynomial")
    private int srpPolynomialDegree = 5;

    @Parameter(valueSet = {"301", "401", "501", "601", "701", "801", "901", "1001"},
            description = "Number of points for the 'flat earth phase' polynomial estimation",
            defaultValue = "501",
            label = "Number of \"Flat Earth\" estimation points")
    private int srpNumberPoints = 501;

    @Parameter(valueSet = {"1", "2", "3", "4", "5"},
            description = "Degree of orbit (polynomial) interpolator",
            defaultValue = "3",
            label = "Orbit interpolation degree")
    private int orbitDegree = 3;

    @Parameter(description = "Use ground square pixel", defaultValue = "true", label = "Square Pixel")
    private Boolean squarePixel = true;

    @Parameter(defaultValue = "false", label = "Subtract topographic phase")
    private boolean subtractTopographicPhase = false;
    /*
        @Parameter(interval = "(1, 10]",
                description = "Degree of orbit interpolation polynomial",
                defaultValue = "3",
                label = "Orbit Interpolation Degree")
        private int orbitDegree = 3;
    */
    @Parameter(description = "The digital elevation model.",
            defaultValue = "SRTM 3Sec",
            label = "Digital Elevation Model")
    private String demName = "SRTM 3Sec";

    @Parameter(label = "External DEM")
    private File externalDEMFile = null;

    @Parameter(label = "DEM No Data Value", defaultValue = "0")
    private double externalDEMNoDataValue = 0;

    @Parameter(label = "External DEM Apply EGM", defaultValue = "true")
    private Boolean externalDEMApplyEGM = true;

    @Parameter(label = "Tile Extension [%]",
            description = "Define extension of tile for DEM simulation (optimization parameter).",
            defaultValue = "100")
    private String tileExtensionPercent = "100";

    // source
    private Map<String, CplxContainer> masterMap = new HashMap<>();
    private Map<String, CplxContainer> slaveMap = new HashMap<>();

    private String[] polarisations;
    private String[] subswaths = new String[]{""};

    // target
    private Map<String, ProductContainer> targetMap = new HashMap<>();
    private Map<Band, Band> detectedSlaveMap = new HashMap<>();

    private boolean isComplex;
    private boolean isTOPSARBurstProduct = false;
    private String productTag = null;
    private Sentinel1Utils su = null;
    private Sentinel1Utils.SubSwathInfo[] subSwath = null;
    private int numSubSwaths = 0;
    private int subSwathIndex = 0;

    private MetadataElement mstRoot = null;
    private MetadataElement slvRoot = null;
    private org.jlinda.core.Point[] mstSceneCentreXYZ = null;
    private HashMap<String, DoubleMatrix> flatEarthPolyMap = new HashMap<>();
    private int sourceImageWidth;
    private int sourceImageHeight;

    private ElevationModel dem = null;
    private double demNoDataValue = 0;
    private double demSamplingLat;
    private double demSamplingLon;

    private static final int ORBIT_DEGREE = 3; // hardcoded
    private static final String PRODUCT_SUFFIX = "_Coh";
    private static final boolean OUTPUT_PHASE = false;
    private static final String FLAT_EARTH_PHASE = "flat_earth_phase";
    private static final String TOPO_PHASE = "topo_phase";


    /**
     * Initializes this operator and sets the one and only target product.
     * <p>The target product can be either defined by a field of type {@link Product} annotated with the
     * {@link TargetProduct TargetProduct} annotation or
     * by calling {@link #setTargetProduct} method.</p>
     * <p>The framework calls this method after it has created this operator.
     * Any client code that must be performed before computation of tile data
     * should be placed here.</p>
     *
     * @throws OperatorException If an error occurs during operator initialisation.
     * @see #getTargetProduct()
     */
    @Override
    public void initialize() throws OperatorException {

        try {
            productTag = "coh";

            mstRoot = AbstractMetadata.getAbstractedMetadata(sourceProduct);
            final MetadataElement slaveElem =
                    sourceProduct.getMetadataRoot().getElement(AbstractMetadata.SLAVE_METADATA_ROOT);
            if (slaveElem != null) {
                slvRoot = slaveElem.getElements()[0];
            }

            checkUserInput();

            constructSourceMetadata();

            constructTargetMetadata();

            createTargetProduct();

            if (isComplex && subtractFlatEarthPhase) {
                if (isTOPSARBurstProduct) {
                    getMstApproxSceneCentreXYZ();
                    constructFlatEarthPolynomialsForTOPSARProduct();
                } else {
                    constructFlatEarthPolynomials();
                }
            }

            if (isComplex && subtractTopographicPhase) {
                defineDEM();
            }
            System.out.println("using box of " + this.cohWinAz + " " + this.cohWinRg);
        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }

    private void checkUserInput() {

        try {
            final InputProductValidator validator = new InputProductValidator(sourceProduct);
            validator.checkIfSARProduct();
            validator.checkIfCoregisteredStack();
            isTOPSARBurstProduct = validator.isTOPSARProduct() && !validator.isDebursted();

            isComplex = AbstractMetadata.getAbstractedMetadata(sourceProduct).
                    getAttributeString(AbstractMetadata.SAMPLE_TYPE).contains("COMPLEX");

            if (isTOPSARBurstProduct) {
                su = new Sentinel1Utils(sourceProduct);
                subswaths = su.getSubSwathNames();
                subSwath = su.getSubSwath();
                numSubSwaths = su.getNumOfSubSwath();
                subSwathIndex = 1; // subSwathIndex is always 1 because of split product
            }

            final String[] polarisationsInBandNames = OperatorUtils.getPolarisations(sourceProduct);
            polarisations = InterferogramOp.getPolsSharedByMstSlv(sourceProduct, polarisationsInBandNames);

            sourceImageWidth = sourceProduct.getSceneRasterWidth();
            sourceImageHeight = sourceProduct.getSceneRasterHeight();
        } catch (Exception e) {
            throw new OperatorException(e);
        }
    }

    private void constructSourceMetadata() throws Exception {

        // define sourceMaster/sourceSlave name tags
        final String masterTag = "mst";
        final String slaveTag = "slv";

        // get sourceMaster & sourceSlave MetadataElement

        // put sourceMaster metadata into the masterMap
        metaMapPut(masterTag, mstRoot, sourceProduct, masterMap);

        // plug sourceSlave metadata into slaveMap
        MetadataElement slaveElem = sourceProduct.getMetadataRoot().getElement(AbstractMetadata.SLAVE_METADATA_ROOT);
        if (slaveElem == null) {
            slaveElem = sourceProduct.getMetadataRoot().getElement("Slave Metadata");
        }
        MetadataElement[] slaveRoot = slaveElem.getElements();
        for (MetadataElement meta : slaveRoot) {
            if (!meta.getName().contains(AbstractMetadata.ORIGINAL_PRODUCT_METADATA)) {
                metaMapPut(slaveTag, meta, sourceProduct, slaveMap);
            }
        }
    }

    private void metaMapPut(final String tag,
                            final MetadataElement root,
                            final Product product,
                            final Map<String, CplxContainer> map) throws Exception {

        for (String swath : subswaths) {
            final String subswath = swath.isEmpty() ? "" : '_' + swath.toUpperCase();

            for (String polarisation : polarisations) {
                final String pol = polarisation.isEmpty() ? "" : '_' + polarisation.toUpperCase();

                // map key: ORBIT NUMBER
                String mapKey = root.getAttributeInt(AbstractMetadata.ABS_ORBIT) + subswath + pol;

                // metadata: construct classes and define bands
                final String date = OperatorUtils.getAcquisitionDate(root);
                final SLCImage meta = new SLCImage(root, product);
                final Orbit orbit = new Orbit(root, ORBIT_DEGREE);
                Band bandReal = null;
                Band bandImag = null;

                // loop through all band names(!) : and pull out only one that matches criteria
                for (String bandName : product.getBandNames()) {
                    if (bandName.contains(tag) && bandName.contains(date)) {
                        if (subswath.isEmpty() || bandName.contains(subswath)) {
                            if (pol.isEmpty() || bandName.contains(pol)) {
                                final Band band = product.getBand(bandName);
                                if (BandUtilsDoris.isBandReal(band)) {
                                    bandReal = band;
                                } else if (BandUtilsDoris.isBandImag(band)) {
                                    bandImag = band;
                                }
                            }
                        }
                    }
                }
                if (bandReal != null && bandImag != null) {
                    map.put(mapKey, new CplxContainer(date, meta, orbit, bandReal, bandImag));
                }
            }
        }
    }

    private void constructTargetMetadata() {

        for (String keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            for (String keySlave : slaveMap.keySet()) {
                final CplxContainer slave = slaveMap.get(keySlave);

                if ((master.polarisation == null || slave.polarisation == null) ||
                        (master.polarisation != null && slave.polarisation != null &&
                                master.polarisation.equals(slave.polarisation))) {
                    // generate name for product bands
                    final String productName = keyMaster + '_' + keySlave;

                    final ProductContainer product = new ProductContainer(productName, master, slave, false);

                    // put ifg-product bands into map
                    targetMap.put(productName, product);
                }
            }
        }
    }

    private void createTargetProduct() {

        targetProduct = new Product(sourceProduct.getName() + PRODUCT_SUFFIX,
                sourceProduct.getProductType(),
                sourceProduct.getSceneRasterWidth(),
                sourceProduct.getSceneRasterHeight());

        ProductUtils.copyProductNodes(sourceProduct, targetProduct);

        if (isComplex) {
            for (String key : targetMap.keySet()) {
                final java.util.List<String> targetBandNames = new ArrayList<>();

                final ProductContainer container = targetMap.get(key);
                final CplxContainer master = container.sourceMaster;
                final CplxContainer slave = container.sourceSlave;

                final String subswath = master.subswath.isEmpty() ? "" : '_' + master.subswath.toUpperCase();
                final String pol = InterferogramOp.getPolarisationTag(master);
                final String tag = subswath + pol + '_' + master.date + '_' + slave.date;

                final String coherenceBandName = productTag + tag;
                final Band coherenceBand = targetProduct.addBand(coherenceBandName, ProductData.TYPE_FLOAT32);
                coherenceBand.setNoDataValueUsed(true);
                coherenceBand.setNoDataValue(master.realBand.getNoDataValue());
                container.addBand(Unit.COHERENCE, coherenceBand.getName());
                coherenceBand.setUnit(Unit.COHERENCE);
                targetBandNames.add(coherenceBand.getName());


                if (subtractTopographicPhase && OUTPUT_PHASE) {
                    final String targetBandTgp = "tgp" + tag;
                    final Band tgpBand = targetProduct.addBand(targetBandTgp, ProductData.TYPE_FLOAT32);
                    container.addBand(Unit.PHASE, tgpBand.getName());
                    tgpBand.setUnit(Unit.PHASE);
                    targetBandNames.add(tgpBand.getName());
                }

                if (subtractFlatEarthPhase && OUTPUT_PHASE) {
                    final String targetBandFep = "fep" + tag;
                    final Band fepBand = targetProduct.addBand(targetBandFep, ProductData.TYPE_FLOAT32);
                    container.addBand(Unit.PHASE, fepBand.getName());
                    fepBand.setUnit(Unit.PHASE);
                    targetBandNames.add(fepBand.getName());
                }

                String slvProductName = StackUtils.findOriginalSlaveProductName(sourceProduct, container.sourceSlave.realBand);
                StackUtils.saveSlaveProductBandNames(targetProduct, slvProductName,
                        targetBandNames.toArray(new String[targetBandNames.size()]));
            }
        } else {
            final int numSrcBands = sourceProduct.getNumBands();
            String[] bandNames = sourceProduct.getBandNames();
            if (numSrcBands < 2) {
                throw new OperatorException("To create a coherence image, more than 2 bands are needed.");
            }
            //masterBand = sourceProduct.getBand(findBandName(bandNames, "mst"));
            //addTargetBand(masterBand.getName(), masterBand.getDataType(), masterBand.getUnit());

            // add slave and coherence bands
            for (int i = 1; i <= numSrcBands; i++) {

                final String slaveBandName = findBandName(bandNames, "slv" + i);
                if (slaveBandName == null) {
                    break;
                }
                final Band slaveBand = sourceProduct.getBand(slaveBandName);
                //addTargetBand(slaveBandName, slaveBand.getDataType(), slaveBand.getUnit());

                final Band coherenceBand = targetProduct.addBand("Coherence_slv" + i, ProductData.TYPE_FLOAT32);
                coherenceBand.setUnit("coherence");
                detectedSlaveMap.put(coherenceBand, slaveBand);
            }
        }
    }

    private static String findBandName(String[] bandNames, String namePattern) {

        String bandName = null;
        for (String name : bandNames) {
            if (name.contains(namePattern)) {
                bandName = name;
                break;
            }
        }
        return bandName;
    }

    private void getMstApproxSceneCentreXYZ() throws Exception {

        final int numOfBursts = subSwath[subSwathIndex - 1].numOfBursts;
        mstSceneCentreXYZ = new org.jlinda.core.Point[numOfBursts];

        for (int b = 0; b < numOfBursts; b++) {
            final double firstLineTime = subSwath[subSwathIndex - 1].burstFirstLineTime[b];
            final double lastLineTime = subSwath[subSwathIndex - 1].burstLastLineTime[b];
            final double slrTimeToFirstPixel = subSwath[subSwathIndex - 1].slrTimeToFirstPixel;
            final double slrTimeToLastPixel = subSwath[subSwathIndex - 1].slrTimeToLastPixel;
            final double latUL = su.getLatitude(firstLineTime, slrTimeToFirstPixel, subSwathIndex);
            final double latUR = su.getLatitude(firstLineTime, slrTimeToLastPixel, subSwathIndex);
            final double latLL = su.getLatitude(lastLineTime, slrTimeToFirstPixel, subSwathIndex);
            final double latLR = su.getLatitude(lastLineTime, slrTimeToLastPixel, subSwathIndex);

            final double lonUL = su.getLongitude(firstLineTime, slrTimeToFirstPixel, subSwathIndex);
            final double lonUR = su.getLongitude(firstLineTime, slrTimeToLastPixel, subSwathIndex);
            final double lonLL = su.getLongitude(lastLineTime, slrTimeToFirstPixel, subSwathIndex);
            final double lonLR = su.getLongitude(lastLineTime, slrTimeToLastPixel, subSwathIndex);

            final double lat = (latUL + latUR + latLL + latLR) / 4.0;
            final double lon = (lonUL + lonUR + lonLL + lonLR) / 4.0;

            final PosVector mstSceneCenter = new PosVector();
            GeoUtils.geo2xyzWGS84(lat, lon, 0.0, mstSceneCenter);
            mstSceneCentreXYZ[b] = new Point(mstSceneCenter.toArray());
        }
    }

    private void constructFlatEarthPolynomialsForTOPSARProduct() throws Exception {

        for (String keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            for (String keySlave : slaveMap.keySet()) {

                CplxContainer slave = slaveMap.get(keySlave);

                for (int s = 0; s < numSubSwaths; s++) {

                    final int numBursts = subSwath[s].numOfBursts;

                    for (int b = 0; b < numBursts; b++) {

                        final String polynomialName = slave.name + '_' + s + '_' + b;

                        flatEarthPolyMap.put(polynomialName, InterferogramOp.estimateFlatEarthPolynomial(
                                master, slave, s + 1, b, mstSceneCentreXYZ, orbitDegree, srpPolynomialDegree,
                                srpNumberPoints, subSwath, su));
                    }
                }
            }
        }
    }

    private void constructFlatEarthPolynomials() throws Exception {

        for (String keyMaster : masterMap.keySet()) {

            CplxContainer master = masterMap.get(keyMaster);

            for (String keySlave : slaveMap.keySet()) {

                CplxContainer slave = slaveMap.get(keySlave);

                flatEarthPolyMap.put(slave.name, InterferogramOp.estimateFlatEarthPolynomial(
                        master.metaData, master.orbit, slave.metaData, slave.orbit, sourceImageWidth,
                        sourceImageHeight, srpPolynomialDegree, srpNumberPoints, sourceProduct));
            }
        }
    }

    private void defineDEM() throws IOException {

        String demResamplingMethod = ResamplingFactory.BILINEAR_INTERPOLATION_NAME;

        if (externalDEMFile == null) {
            dem = DEMFactory.createElevationModel(demName, demResamplingMethod);
            demNoDataValue = dem.getDescriptor().getNoDataValue();
            demSamplingLat = dem.getDescriptor().getTileWidthInDegrees() * (1.0f /
                    dem.getDescriptor().getTileWidth()) * org.jlinda.core.Constants.DTOR;

            demSamplingLon = demSamplingLat;

        } else {

            dem = new FileElevationModel(externalDEMFile, demResamplingMethod, externalDEMNoDataValue);
            ((FileElevationModel) dem).applyEarthGravitionalModel(externalDEMApplyEGM);
            demNoDataValue = externalDEMNoDataValue;
            demName = externalDEMFile.getName();

            try {
                demSamplingLat =
                        (dem.getGeoPos(new PixelPos(0, 1)).getLat() - dem.getGeoPos(new PixelPos(0, 0)).getLat()) *
                                org.jlinda.core.Constants.DTOR;
                demSamplingLon =
                        (dem.getGeoPos(new PixelPos(1, 0)).getLon() - dem.getGeoPos(new PixelPos(0, 0)).getLon()) *
                                org.jlinda.core.Constants.DTOR;
            } catch (Exception e) {
                throw new OperatorException("The DEM '" + demName + "' cannot be properly interpreted.");
            }
        }
    }

    @Override
    public void computeTileStack(Map<Band, Tile> targetTileMap, Rectangle targetRectangle, ProgressMonitor pm)
            throws OperatorException {

        if (isTOPSARBurstProduct) {
            computeTileForTOPSARProduct(targetTileMap, targetRectangle, pm);
        }
    }

    private void computeTileForTOPSARProduct(
            final Map<Band, Tile> targetTileMap, final Rectangle targetRectangle, final ProgressMonitor pm)
            throws OperatorException {

        try {
            final int tx0 = targetRectangle.x;
            final int ty0 = targetRectangle.y;
            final int tw = targetRectangle.width;
            final int th = targetRectangle.height;
            final int txMax = tx0 + tw;
            final int tyMax = ty0 + th;
            //System.out.println("tx0 = " + tx0 + ", ty0 = " + ty0 + ", tw = " + tw + ", th = " + th);

            for (int burstIndex = 0; burstIndex < subSwath[subSwathIndex - 1].numOfBursts; burstIndex++) {
                final int firstLineIdx = burstIndex * subSwath[subSwathIndex - 1].linesPerBurst;
                final int lastLineIdx = firstLineIdx + subSwath[subSwathIndex - 1].linesPerBurst - 1;

                if (tyMax <= firstLineIdx || ty0 > lastLineIdx) {
                    continue;
                }

                final int ntx0 = tx0;
                final int ntw = tw;
                final int nty0 = Math.max(ty0, firstLineIdx);
                final int ntyMax = Math.min(tyMax, lastLineIdx + 1);
                final int nth = ntyMax - nty0;
                final Rectangle partialTileRectangle = new Rectangle(ntx0, nty0, ntw, nth);
                //System.out.println("burst = " + burstIndex + ": ntx0 = " + ntx0 + ", nty0 = " + nty0 + ", ntw = " + ntw + ", nth = " + nth);

                computePartialTile(subSwathIndex, burstIndex, firstLineIdx, partialTileRectangle, targetTileMap);
            }

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        } finally {
            pm.done();
        }
    }

    private void computePartialTile(final int subSwathIndex, final int burstIndex,
                                    final int firstLineIdx, final Rectangle targetRectangle,
                                    final Map<Band, Tile> targetTileMap) throws OperatorException {

        try {
            final BorderExtender border = BorderExtender.createInstance(BorderExtender.BORDER_ZERO);

            final int y0 = targetRectangle.y;
            final int yN = y0 + targetRectangle.height - 1;
            final int x0 = targetRectangle.x;
            final int xN = x0 + targetRectangle.width - 1;

            final int cohx0 = targetRectangle.x - (cohWinRg - 1) / 2;
            final int cohy0 = targetRectangle.y - (cohWinAz - 1) / 2;
            final int cohw = targetRectangle.width + cohWinRg - 1;
            final int cohh = targetRectangle.height + cohWinAz - 1;
            final Rectangle extRect = new Rectangle(cohx0, cohy0, cohw, cohh);

            final org.jlinda.core.Window tileWindow = new org.jlinda.core.Window(
                    cohy0 - firstLineIdx, cohy0 + cohh - 1 - firstLineIdx, cohx0, cohx0 + cohw - 1);

            final SLCImage mstMeta = targetMap.values().iterator().next().sourceMaster.metaData.clone();
            updateMstMetaData(burstIndex, mstMeta);
            final Orbit mstOrbit = targetMap.values().iterator().next().sourceMaster.orbit;

            DemTile demTile = null;


            final int minLine = 0;
            final int maxLine = subSwath[subSwathIndex - 1].linesPerBurst - 1;
            final int minPixel = 0;
            final int maxPixel = subSwath[subSwathIndex - 1].samplesPerBurst - 1;

            for (String cohKey : targetMap.keySet()) {

                final ProductContainer product = targetMap.get(cohKey);
                final SLCImage slvMeta = product.sourceSlave.metaData.clone();
                updateSlvMetaData(product, burstIndex, slvMeta);
                final Orbit slvOrbit = product.sourceSlave.orbit;

                final Tile mstTileReal = getSourceTile(product.sourceMaster.realBand, extRect, border);
                final Tile mstTileImag = getSourceTile(product.sourceMaster.imagBand, extRect, border);
                final ComplexDoubleMatrix dataMaster = TileUtilsDoris.pullComplexDoubleMatrix(mstTileReal, mstTileImag);

                final Tile slvTileReal = getSourceTile(product.sourceSlave.realBand, extRect, border);
                final Tile slvTileImag = getSourceTile(product.sourceSlave.imagBand, extRect, border);
                final ComplexDoubleMatrix dataSlave = TileUtilsDoris.pullComplexDoubleMatrix(slvTileReal, slvTileImag);

                final String polynomialName = product.sourceSlave.name + '_' + (subSwathIndex - 1) + '_' + burstIndex;


                for (int i = 0; i < dataMaster.length; i++) {
                    double tmp = norm(dataMaster.get(i));
                    dataMaster.put(i, dataMaster.get(i).mul(dataSlave.get(i).conj()));
                    dataSlave.put(i, new ComplexDouble(norm(dataSlave.get(i)), tmp));
                }

                DoubleMatrix cohMatrix = coherence2(dataMaster, dataSlave, cohWinAz, cohWinRg);

                saveCoherence(cohMatrix, product, targetTileMap, targetRectangle);
            }

        } catch (Throwable e) {
            OperatorUtils.catchOperatorException(getId(), e);
        }
    }


    private double norm(final ComplexDouble number) {
        return number.real() * number.real() + number.imag() * number.imag();
    }

    public DoubleMatrix coherence2(final ComplexDoubleMatrix input, final ComplexDoubleMatrix norms, final int winL, final int winP) {


        final int extent_RG = input.columns;
        final int extent_AZ = input.rows - winL + 1;
        final DoubleMatrix result = new DoubleMatrix(input.rows - winL + 1, input.columns - winP + 1);

        // temp variables
        int i, j, k, l;
        ComplexDouble sum;
        ComplexDouble power;
        final int leadingZeros = (winP - 1) / 2;  // number of pixels=0 floor...
        final int trailingZeros = (winP) / 2;     // floor...

        for (j = leadingZeros; j < extent_RG - trailingZeros; j++) {

            sum = new ComplexDouble(0);
            power = new ComplexDouble(0);

            //// Compute sum over first data block ////
            int minL = j - leadingZeros;
            int maxL = minL + winP;
            for (k = 0; k < winL; k++) {
                for (l = minL; l < maxL; l++) {
                    //sum.addi(input.get(k, l));
                    //power.addi(norms.get(k, l));
                    int inI = 2 * input.index(k, l);
                    sum.set(sum.real() + input.data[inI], sum.imag() + input.data[inI + 1]);
                    power.set(power.real() + norms.data[inI], power.imag() + norms.data[inI + 1]);
                }
            }
            result.put(0, minL, coherenceProduct(sum, power));

            //// Compute (relatively) sum over rest of data blocks ////
            final int maxI = extent_AZ - 1;
            for (i = 0; i < maxI; i++) {
                final int iwinL = i + winL;
                for (l = minL; l < maxL; l++) {
                    //sum.addi(input.get(iwinL, l).sub(input.get(i, l)));
                    //power.addi(norms.get(iwinL, l).sub(norms.get(i, l)));

                    int inI = 2 * input.index(i, l);
                    int inWinL = 2 * input.index(iwinL, l);
                    sum.set(sum.real() + (input.data[inWinL] - input.data[inI]), sum.imag() + (input.data[inWinL + 1] - input.data[inI + 1]));
                    power.set(power.real() + (norms.data[inWinL] - norms.data[inI]), power.imag() + (norms.data[inWinL + 1] - norms.data[inI + 1]));
                }
                result.put(i + 1, j - leadingZeros, coherenceProduct(sum, power));
            }
        }
        return result;
    }

    static double coherenceProduct(final ComplexDouble sum, final ComplexDouble power) {
        final double product = power.real() * power.imag();
        return (product > 0.0) ? sum.abs() / Math.sqrt(product) : 0.0;
    }

    private void saveCoherence(final DoubleMatrix cohMatrix, final ProductContainer product,
                               final Map<Band, Tile> targetTileMap, final Rectangle targetRectangle) {

        final int x0 = targetRectangle.x;
        final int y0 = targetRectangle.y;
        final int maxX = x0 + targetRectangle.width;
        final int maxY = y0 + targetRectangle.height;

        final Band coherenceBand = targetProduct.getBand(product.getBandName(Unit.COHERENCE));
        final Tile coherenceTile = targetTileMap.get(coherenceBand);
        final ProductData coherenceData = coherenceTile.getDataBuffer();

        final double srcNoDataValue = product.sourceMaster.realBand.getNoDataValue();
        final Tile slvTileReal = getSourceTile(product.sourceSlave.realBand, targetRectangle);
        final ProductData srcSlvData = slvTileReal.getDataBuffer();
        final TileIndex srcSlvIndex = new TileIndex(slvTileReal);

        final TileIndex tgtIndex = new TileIndex(coherenceTile);
        for (int y = y0; y < maxY; y++) {
            tgtIndex.calculateStride(y);
            srcSlvIndex.calculateStride(y);
            final int yy = y - y0;
            for (int x = x0; x < maxX; x++) {
                final int tgtIdx = tgtIndex.getIndex(x);
                final int xx = x - x0;

                if (srcSlvData.getElemDoubleAt(srcSlvIndex.getIndex(x)) == srcNoDataValue) {
                    coherenceData.setElemFloatAt(tgtIdx, (float) srcNoDataValue);
                } else {
                    final double coh = cohMatrix.get(yy, xx);
                    coherenceData.setElemFloatAt(tgtIdx, (float) coh);
                }
            }
        }
    }


    //-----------------------------------------------------------------------------
    //
    //                          META DATA UPDATES
    //
    //-----------------------------------------------------------------------------
    private void updateMstMetaData(final int burstIndex, final SLCImage mstMeta) {

        final double burstFirstLineTimeMJD = subSwath[subSwathIndex - 1].burstFirstLineTime[burstIndex] /
                Constants.secondsInDay;

        final double burstFirstLineTimeSecondsOfDay = (burstFirstLineTimeMJD - (int)burstFirstLineTimeMJD) *
                Constants.secondsInDay;

        mstMeta.settAzi1(burstFirstLineTimeSecondsOfDay);

        mstMeta.setCurrentWindow(new org.jlinda.core.Window(0, subSwath[subSwathIndex - 1].linesPerBurst - 1,
                0, subSwath[subSwathIndex - 1].samplesPerBurst - 1));

        mstMeta.setOriginalWindow(new org.jlinda.core.Window(0, subSwath[subSwathIndex - 1].linesPerBurst - 1,
                0, subSwath[subSwathIndex - 1].samplesPerBurst - 1));

        mstMeta.setApproxGeoCentreOriginal(getApproxGeoCentre(subSwathIndex, burstIndex));
    }

    private void updateSlvMetaData(final ProductContainer product, final int burstIndex, final SLCImage slvMeta) {

        final double slvBurstFirstLineTimeMJD = slvMeta.getMjd() - product.sourceMaster.metaData.getMjd() +
                subSwath[subSwathIndex - 1].burstFirstLineTime[burstIndex] / Constants.secondsInDay;

        final double slvBurstFirstLineTimeSecondsOfDay = (slvBurstFirstLineTimeMJD - (int)slvBurstFirstLineTimeMJD) *
                Constants.secondsInDay;

        slvMeta.settAzi1(slvBurstFirstLineTimeSecondsOfDay);

        slvMeta.setCurrentWindow(new org.jlinda.core.Window(0, subSwath[subSwathIndex - 1].linesPerBurst - 1,
                0, subSwath[subSwathIndex - 1].samplesPerBurst - 1));

        slvMeta.setOriginalWindow(new org.jlinda.core.Window(0, subSwath[subSwathIndex - 1].linesPerBurst - 1,
                0, subSwath[subSwathIndex - 1].samplesPerBurst - 1));
    }

    private GeoPoint getApproxGeoCentre(final int subSwathIndex, final int burstIndex) {

        final int cols = subSwath[subSwathIndex - 1].latitude[0].length;

        double lat = 0.0, lon = 0.0;
        for (int r = burstIndex; r <= burstIndex + 1; r++) {
            for (int c = 0; c < cols; c++) {
                lat += subSwath[subSwathIndex - 1].latitude[r][c];
                lon += subSwath[subSwathIndex - 1].longitude[r][c];
            }
        }

        return new GeoPoint(lat / (2*cols), lon / (2*cols));
    }
}
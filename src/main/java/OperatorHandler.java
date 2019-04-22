import com.bc.ceres.core.ProgressMonitor;
import operators.CoherenceOp;
import operators.TestCoherenceOp;
import org.esa.snap.core.dataio.ProductIOPlugInManager;
import org.esa.snap.core.dataio.dimap.DimapProductReader;
import org.esa.snap.core.dataio.dimap.DimapProductReaderPlugIn;
import org.esa.snap.core.dataio.dimap.DimapProductWriter;
import org.esa.snap.core.dataio.dimap.DimapProductWriterPlugIn;
import org.esa.snap.core.datamodel.Band;
import org.esa.snap.core.datamodel.Product;
import org.esa.snap.core.gpf.Operator;
import org.esa.snap.core.gpf.Tile;
import org.esa.snap.core.gpf.internal.TileImpl;
import org.esa.snap.core.image.ImageManager;
import org.esa.snap.core.util.ImageUtils;
import org.esa.snap.dataio.geotiff.GeoTiffProductReaderPlugIn;
import org.esa.snap.dataio.geotiff.GeoTiffProductWriterPlugIn;

import javax.media.jai.RasterFactory;
import java.awt.*;
import java.awt.image.SampleModel;
import java.awt.image.WritableRaster;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

public class OperatorHandler {

    private String outputDirectory;
    private String outputProductName;
    protected DimapProductReaderPlugIn readerPlugIn1;
    protected DimapProductWriterPlugIn writerPlugIn1;

    protected DimapProductReader reader1;

    protected DimapProductWriter writer1;

    protected int productWidth;
    protected int productHeight;
    protected int tileWidth;
    protected int tileHeight;

    public OperatorHandler(String outputDirectory, String outputProductName) {
        this.outputDirectory = outputDirectory;
        this.outputProductName = outputProductName;
        this.readerPlugIn1 = new DimapProductReaderPlugIn();
        this.writerPlugIn1 = new DimapProductWriterPlugIn();
    }

    private void initio(){
        this.reader1 = new DimapProductReader(this.readerPlugIn1);
        this.writer1 = new DimapProductWriter(this.writerPlugIn1);
        ProductIOPlugInManager.getInstance().addReaderPlugIn(new DimapProductReaderPlugIn());
        ProductIOPlugInManager.getInstance().addReaderPlugIn(new GeoTiffProductReaderPlugIn());
        ProductIOPlugInManager.getInstance().addWriterPlugIn(new DimapProductWriterPlugIn());
        ProductIOPlugInManager.getInstance().addWriterPlugIn(new GeoTiffProductWriterPlugIn());
    }

    private void closeio() {
        try {
            this.reader1.close();

        } catch (IOException e) {
            System.out.println("Can not close reader 1");
            e.printStackTrace();
        }

        try {
            this.writer1.close();
        } catch (IOException e) {
            System.out.println("Can not close writer 1");
            e.printStackTrace();
        }
    }

    public void coherenceEstimation(String product1Path){
        this.initio();
        File file1 = new File(product1Path);

        File outputDir = new File(this.outputDirectory);
        File outputProduct = new File(this.outputProductName);
        try {
            Product product1 = this.reader1.readProductNodes(file1, null);
            this.tileHeight = 560;
            this.tileWidth = 404;

            Operator operator = new CoherenceOp();
            operator.setSourceProduct(product1);
            for(String name: product1.getBandNames()){
                System.out.println(name);
            }
            System.out.println(product1.getNumBands());
            //((CoherenceOp) operator).sourceProduct = product1;
            operator.initialize();
            Product target = operator.getTargetProduct();
            target.setPreferredTileSize(this.tileWidth,this.tileHeight);
            System.out.println("----------------target product---------------");
            System.out.println(target.getName());
            System.out.println(target.getNumBands());
            for(String name :target.getBandNames()){
                System.out.println(name);
            }


            this.writer1.initDirs(outputDir);
            this.writer1.writeProductNodes(target,outputProduct);
            target.setPreferredTileSize(this.tileWidth,this.tileHeight);
            Band cohBand = target.getBand("coh_IW1_VV_02Jun2017_14Jun2017");

            this.productHeight = cohBand.getRasterHeight();
            this.productWidth = cohBand.getRasterWidth();

            long begin = System.currentTimeMillis();
            int currentX=0, currentY=0;
            /*for(currentY =0; currentY < this.productHeight; currentY += this.tileHeight) {
                for (currentX = 0; currentX < this.productWidth; currentX += this.tileWidth) {*/
                    List<Thread> workers = new ArrayList<>();
                    Rectangle currentRectangle = this.getCheckedTile(currentX, currentY, this.tileWidth, this.tileHeight);
                    int cohDataType = ImageManager.getDataBufferType(cohBand.getDataType());

                    SampleModel cohSample = ImageUtils.createSingleBandedSampleModel(cohDataType, currentRectangle.width, currentRectangle.height);
                    WritableRaster cohRaster = RasterFactory.createWritableRaster(cohSample, new Point(currentX, currentY));

                    Map<Band, Tile> tiles = new HashMap<>();
                    Tile cohTile = new TileImpl(cohBand, cohRaster);
                    tiles.put(cohBand, cohTile);

                    operator.computeTileStack(tiles, currentRectangle, ProgressMonitor.NULL);

                    tiles.forEach((band, tile) -> {
                        Thread tempThread = new WriterThread(band,currentRectangle,tile.getDataBuffer(), this.writer1);
                        tempThread.start();
                        workers.add(tempThread);
                    });

                    /*Thread tempThread = new CopyThread(isource,itarget,currentRectangle,this.reader1,this.writer1);
                    tempThread.start();
                    workers.add(tempThread);
                    tempThread = new CopyThread(qsource,qtarget,currentRectangle,this.reader1,this.writer1);
                    tempThread.start();
                    workers.add(tempThread);*/

                    for (Thread worker : workers) {
                        worker.join();
                    }
            /*    }
            }*/
            System.out.println("Raster calc took: " + (System.currentTimeMillis() - begin));
            /*System.out.println("slave part took: " + GPUBackGeocodingOp.slaveTotal);
            System.out.println("deramp demod part took: " + GPUBackGeocodingOp.derampTotal);
            System.out.println("interpolation part took: " + GPUBackGeocodingOp.interpolationTotal);*/
        }catch (Exception ioe){
            System.out.println("problems reading nodes");
            ioe.printStackTrace();
        }

        this.closeio();
    }

    protected Rectangle getCheckedTile(int offsetX, int offsetY, int width, int height){
        if((offsetY + this.tileHeight) > this.productHeight){
            height = this.productHeight - offsetY;
        }
        if((offsetX + this.tileWidth) > this.productWidth){
            width = this.productWidth - offsetX;
        }
        return new Rectangle(offsetX, offsetY, width, height);
    }
}

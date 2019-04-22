import com.bc.ceres.core.ProgressMonitor;
import org.esa.snap.core.dataio.ProductWriter;
import org.esa.snap.core.datamodel.Band;
import org.esa.snap.core.datamodel.ProductData;

import java.awt.*;
import java.io.IOException;

public class WriterThread extends Thread {

    private Band band;
    private Rectangle rectangle;
    private ProductData data;
    private ProductWriter writer;

    public WriterThread(Band band, Rectangle rectangle, ProductData data, ProductWriter writer) {
        super();
        this.band = band;
        this.rectangle = rectangle;
        this.data = data;
        this.writer = writer;
    }

    @Override
    public void run() {

        try {
            this.writer.writeBandRasterData(this.band,
                    this.rectangle.x,
                    this.rectangle.y,
                    this.rectangle.width,
                    this.rectangle.height,
                    this.data,
                    ProgressMonitor.NULL);

        } catch (IOException e) {
            System.out.println("Can not write for band " + this.band.getName());
            e.printStackTrace();
        }
    }
}

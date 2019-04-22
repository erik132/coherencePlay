public class Main {

    public static void main(String[] args){
        System.out.println("It's Coherence time children!!!");
        String outputDirectory = "C:\\Users\\erik.soekov\\esaSnap\\snapDebugsTests\\coherencePlay";
        String outputProduct = "C:\\Users\\erik.soekov\\esaSnap\\snapDebugsTests\\coherencePlay.dim";
        String inputPath = "C:\\Users\\erik.soekov\\esaSnap\\snapDebugsTests\\malbolgeEtalon.dim";
        OperatorHandler handler = new OperatorHandler(outputDirectory, outputProduct);
        handler.coherenceEstimation(inputPath);

    }
}

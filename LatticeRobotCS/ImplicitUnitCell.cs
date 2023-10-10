using g4;

public class ImplicitUnitCell : BoundedImplicitFunction3d {
    public int LatticeIndex;
    public Dictionary<string, double> Parameters { get; set; }
    
    private MyImplicit myImplicit;

    public ImplicitUnitCell(string code, int latticeIndex) {
        Parameters = new Dictionary<string, double>();
        Parameters["size_x"] = 10.0;
        Parameters["size_y"] = 10.0;
        Parameters["size_z"] = 10.0;
        Parameters["bias"] = 0.0;
        Parameters["thickness"] = 1.0;
        Parameters["drop_x"] = 1.0;
        Parameters["drop_y"] = 1.0;
        Parameters["drop_z"] = 1.0;
        Parameters["gyroid"] = 0.0;

        LatticeIndex = latticeIndex;
        myImplicit = new MyImplicit(latticeIndex, Parameters);
    }

    public double Value(ref Vector3d p) {
        return myImplicit.IndexedLattice(p).Distance;
    }

    public AxisAlignedBox3d Bounds()
    {
        var halfsize = new Vector3d(Parameters["size_x"], Parameters["size_y"], Parameters["size_z"]) * 0.5;
        return new AxisAlignedBox3d(-halfsize, halfsize);
    }


}
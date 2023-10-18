using g4;

public class ImplicitUnitCell : BoundedImplicitFunction3d {
    public int LatticeIndex;
    public Dictionary<string, double> Parameters { get; set; }
    
    private DiamondImplicit diamondImplicit;

    public ImplicitUnitCell(int latticeIndex) {
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
        diamondImplicit = new DiamondImplicit(latticeIndex, Parameters);
    }

    public double Value(ref Vector3d p) {
        return diamondImplicit.IndexedLattice(p).Distance;
    }

    public AxisAlignedBox3d Bounds() {
        var halfsize = new Vector3d(Parameters["size_x"], Parameters["size_y"], Parameters["size_z"]) * 0.5;
        return new AxisAlignedBox3d(-halfsize, halfsize);
    }


}
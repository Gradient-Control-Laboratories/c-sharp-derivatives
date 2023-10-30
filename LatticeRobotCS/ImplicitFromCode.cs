using g4;

public class ImplicitFromCode : BoundedImplicitFunction3d{
    public string Code;
    public AxisAlignedBox3d Box;

    public ImplicitFromCode(string code, AxisAlignedBox3d box) {
        Code = code;
        Box = box;
    }

    public double Value(ref Vector3d p) {
        return Math.Sin(p.x) * Math.Cos(p.y) +
               Math.Sin(p.y) * Math.Cos(p.z) +
               Math.Sin(p.z) * Math.Cos(p.x);
    }

    public AxisAlignedBox3d Bounds()
    {
        return Box;
    }


}
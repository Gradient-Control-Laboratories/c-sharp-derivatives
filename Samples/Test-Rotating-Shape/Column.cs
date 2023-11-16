using g4;
using System.Collections.Generic;

public partial class LRImplicit
{
    public static int VariantIndex;
    public static double size_x;
    public static double size_y;
    public static double size_z;
    public static double bias;
    public static double thickness;
    public static double length;

    static Implicit vdot(Implicit ax, Implicit ay, Implicit az, Implicit bx, Implicit by, Implicit bz)
    {
        Implicit _vdot_000 = Multiply(ax, bx);
        Implicit _vdot_001 = Multiply(ay, by);
        Implicit _vdot_002 = Add(_vdot_000, _vdot_001);
        Implicit _vdot_003 = Multiply(az, bz);
        Implicit _vdot_004 = Add(_vdot_002, _vdot_003);
        return _vdot_004;
    }

    static Implicit vdot(Implicit ax, Implicit ay, Implicit az, double bx, double by, double bz)
    {
        Implicit _vdot_000 = Multiply(ax, bx);
        Implicit _vdot_001 = Multiply(ay, by);
        Implicit _vdot_002 = Add(_vdot_000, _vdot_001);
        Implicit _vdot_003 = Multiply(az, bz);
        Implicit _vdot_004 = Add(_vdot_002, _vdot_003);
        return _vdot_004;
    }

    static double vdot(double ax, double ay, double az, double bx, double by, double bz)
    {
        double _vdot_000 = ax * bx;
        double _vdot_001 = ay * by;
        double _vdot_002 = _vdot_000 + _vdot_001;
        double _vdot_003 = az * bz;
        double _vdot_004 = _vdot_002 + _vdot_003;
        return _vdot_004;
    }

    static Implicit vlength(Implicit dx, Implicit dy, Implicit dz)
    {
        Implicit _vlength_000 = vdot(dx, dy, dz, dx, dy, dz);
        Implicit _vlength_001 = Sqrt(_vlength_000);
        return _vlength_001;
    }

    static Implicit fclamp(Implicit t, double a, double b)
    {
        Implicit _fclamp_000 = Min(t, b);
        Implicit _fclamp_001 = Max(a, _fclamp_000);
        return _fclamp_001;
    }

    static Implicit LineSegment(Implicit px, Implicit py, Implicit pz, double ax, double ay, double az, double bx, double by, double bz)
    {
        Implicit pax = Subtract(px, ax);
        double _bax_000 = bx - ax;
        double bax = _bax_000 + 0.0;
        Implicit pay = Subtract(py, ay);
        Implicit paz = Subtract(pz, az);
        double _bay_000 = by - ay;
        double bay = _bay_000 + 0.0;
        double _baz_000 = bz - az;
        double baz = _baz_000 + 0.0;
        Implicit _h_000 = vdot(pax, pay, paz, bax, bay, baz);
        double _h_001 = vdot(bax, bay, baz, bax, bay, baz);
        Implicit _h_002 = Divide(_h_000, _h_001);
        Implicit h = fclamp(_h_002, 0.0, 1.0);
        Implicit _dx_000 = Multiply(bax, h);
        Implicit dx = Subtract(pax, _dx_000);
        Implicit _dy_000 = Multiply(bay, h);
        Implicit dy = Subtract(pay, _dy_000);
        Implicit _dz_000 = Multiply(baz, h);
        Implicit dz = Subtract(paz, _dz_000);
        Implicit _LineSegment_000 = vlength(dx, dy, dz);
        return _LineSegment_000;
    }

    static Implicit baseLattice(Vector3d p)
    {
        Implicit x = new Implicit(p.x, new Vector3d(1.0, 0.0, 0.0));
        Implicit y = new Implicit(p.y, new Vector3d(0.0, 1.0, 0.0));
        Implicit z = new Implicit(p.z, new Vector3d(0.0, 0.0, 1.0));
        double halfLength = length * 0.5;
        double _beam_000 = -1.0 * halfLength;
        Implicit beam = LineSegment(x, y, z, 0.0, 0.0, _beam_000, 0.0, 0.0, halfLength);
        Implicit lattice = Subtract(beam, bias);
        return lattice;
    }

    public static Implicit IndexedLattice(Vector3d p)
    {
        Implicit lattice = baseLattice(p);
        Implicit solid = lattice;
        if (VariantIndex == 0) return solid;
        Implicit inverse = Multiply(-1.0, lattice);
        if (VariantIndex == 1) return inverse;
        Implicit _thin_000 = Abs(lattice);
        double _thin_001 = thickness * 0.5;
        Implicit thin = Subtract(_thin_000, _thin_001);
        if (VariantIndex == 2) return thin;
        Implicit twin = Multiply(-1.0, thin);
        if (VariantIndex == 3) return twin;
        Implicit unknown = Sphere(p, new Vector3d(0.0), 0.5);
        return unknown;
    }

    public static Implicit ScaledLattice(Vector3d scaledP)
    {
        Vector3d p = (scaledP - center) * 10.0;
        Implicit result = IndexedLattice(p);
        Implicit indexed = Divide(result, 10.0);
        return indexed;
    }

}
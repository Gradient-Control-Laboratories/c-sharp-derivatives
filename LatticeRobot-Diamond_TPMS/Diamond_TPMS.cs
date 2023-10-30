using g4;
using System.Collections.Generic;

public partial class LRImplicit  {
    public static int LatticeIndex;

    public static double size_x;
    public static double size_y;
    public static double size_z;
    public static double bias;
    public static double thickness;
    public static double drop_x;
    public static double drop_y;
    public static double drop_z;
    public static double gyroid;

    public LRImplicit() {}
    public static void Initialize(int latticeIndex, Dictionary<string, ImplicitParameter> parameters) {
        LatticeIndex = latticeIndex;

        size_x = parameters["size_x"].defaultValue;
        size_y = parameters["size_y"].defaultValue;
        size_z = parameters["size_z"].defaultValue;
        bias = parameters["bias"].defaultValue;
        thickness = parameters["thickness"].defaultValue;
        drop_x = parameters["drop_x"].defaultValue;
        drop_y = parameters["drop_y"].defaultValue;
        drop_z = parameters["drop_z"].defaultValue;
        gyroid = parameters["gyroid"].defaultValue;
    }

    public static double Value(ref Vector3d p) {
        return IndexedLattice(p).Distance;
    }

    static double supremum = SQRT2;

    static Implicit mix(Implicit a, double b, double t) {
        double _mix_000 = 1.0 - t;
        Implicit _mix_001 = Multiply(_mix_000, a);
        double _mix_002 = t * b;
        Implicit _mix_003 = Add(_mix_001, _mix_002);
        return _mix_003;
    }

    static Implicit solidLattice(Vector3d p) {
        double _tpms_000 = 1.0 - gyroid;
        Implicit x = new Implicit(p.x, new Vector3d(1.0, 0.0, 0.0));
        double halfsizeX = size_x * 0.5;
        double frequencyX = PI / halfsizeX;
        Implicit _sx_000 = Multiply(x, frequencyX);
        Implicit sx = Sin(_sx_000);
        Implicit _tpms_001 = Multiply(_tpms_000, sx);
        Implicit y = new Implicit(p.y, new Vector3d(0.0, 1.0, 0.0));
        double halfsizeY = size_y * 0.5;
        double frequencyY = PI / halfsizeY;
        Implicit _sy_000 = Multiply(y, frequencyY);
        Implicit sy = Sin(_sy_000);
        Implicit _tpms_002 = Multiply(_tpms_001, sy);
        Implicit z = new Implicit(p.z, new Vector3d(0.0, 0.0, 1.0));
        double halfsizeZ = size_z * 0.5;
        double frequencyZ = PI / halfsizeZ;
        Implicit _sz_000 = Multiply(z, frequencyZ);
        Implicit sz = Sin(_sz_000);
        Implicit _tpms_003 = Multiply(_tpms_002, sz);
        Implicit _tpms_004 = Multiply(drop_x, sx);
        Implicit _cy_000 = Multiply(y, frequencyY);
        Implicit cy = Cos(_cy_000);
        Implicit _tpms_005 = Multiply(_tpms_004, cy);
        Implicit _cz_000 = Multiply(z, frequencyZ);
        Implicit cz = Cos(_cz_000);
        Implicit _tpms_006 = mix(cz, 1.0, gyroid);
        Implicit _tpms_007 = Multiply(_tpms_005, _tpms_006);
        Implicit _tpms_008 = Add(_tpms_003, _tpms_007);
        Implicit _tpms_009 = Multiply(drop_y, sy);
        Implicit _tpms_010 = Multiply(_tpms_009, cz);
        Implicit _cx_000 = Multiply(x, frequencyX);
        Implicit cx = Cos(_cx_000);
        Implicit _tpms_011 = mix(cx, 1.0, gyroid);
        Implicit _tpms_012 = Multiply(_tpms_010, _tpms_011);
        Implicit _tpms_013 = Add(_tpms_008, _tpms_012);
        Implicit _tpms_014 = Multiply(drop_z, sz);
        Implicit _tpms_015 = Multiply(_tpms_014, cx);
        Implicit _tpms_016 = mix(cy, 1.0, gyroid);
        Implicit _tpms_017 = Multiply(_tpms_015, _tpms_016);
        Implicit tpms = Add(_tpms_013, _tpms_017);
        Implicit _lattice_000 = Multiply(tpms, size_x);
        double infimum = -1.0 * SQRT2;
        double range = supremum - infimum;
        double _lattice_001 = 2.0 * range;
        Implicit _lattice_002 = Divide(_lattice_000, _lattice_001);
        Implicit lattice = Subtract(_lattice_002, bias);
        Implicit solid = lattice;
        return solid;
    }

    static Implicit inverseLattice(Vector3d p) {
        double _tpms_000 = 1.0 - gyroid;
        Implicit x = new Implicit(p.x, new Vector3d(1.0, 0.0, 0.0));
        double halfsizeX = size_x * 0.5;
        double frequencyX = PI / halfsizeX;
        Implicit _sx_000 = Multiply(x, frequencyX);
        Implicit sx = Sin(_sx_000);
        Implicit _tpms_001 = Multiply(_tpms_000, sx);
        Implicit y = new Implicit(p.y, new Vector3d(0.0, 1.0, 0.0));
        double halfsizeY = size_y * 0.5;
        double frequencyY = PI / halfsizeY;
        Implicit _sy_000 = Multiply(y, frequencyY);
        Implicit sy = Sin(_sy_000);
        Implicit _tpms_002 = Multiply(_tpms_001, sy);
        Implicit z = new Implicit(p.z, new Vector3d(0.0, 0.0, 1.0));
        double halfsizeZ = size_z * 0.5;
        double frequencyZ = PI / halfsizeZ;
        Implicit _sz_000 = Multiply(z, frequencyZ);
        Implicit sz = Sin(_sz_000);
        Implicit _tpms_003 = Multiply(_tpms_002, sz);
        Implicit _tpms_004 = Multiply(drop_x, sx);
        Implicit _cy_000 = Multiply(y, frequencyY);
        Implicit cy = Cos(_cy_000);
        Implicit _tpms_005 = Multiply(_tpms_004, cy);
        Implicit _cz_000 = Multiply(z, frequencyZ);
        Implicit cz = Cos(_cz_000);
        Implicit _tpms_006 = mix(cz, 1.0, gyroid);
        Implicit _tpms_007 = Multiply(_tpms_005, _tpms_006);
        Implicit _tpms_008 = Add(_tpms_003, _tpms_007);
        Implicit _tpms_009 = Multiply(drop_y, sy);
        Implicit _tpms_010 = Multiply(_tpms_009, cz);
        Implicit _cx_000 = Multiply(x, frequencyX);
        Implicit cx = Cos(_cx_000);
        Implicit _tpms_011 = mix(cx, 1.0, gyroid);
        Implicit _tpms_012 = Multiply(_tpms_010, _tpms_011);
        Implicit _tpms_013 = Add(_tpms_008, _tpms_012);
        Implicit _tpms_014 = Multiply(drop_z, sz);
        Implicit _tpms_015 = Multiply(_tpms_014, cx);
        Implicit _tpms_016 = mix(cy, 1.0, gyroid);
        Implicit _tpms_017 = Multiply(_tpms_015, _tpms_016);
        Implicit tpms = Add(_tpms_013, _tpms_017);
        Implicit _lattice_000 = Multiply(tpms, size_x);
        double infimum = -1.0 * SQRT2;
        double range = supremum - infimum;
        double _lattice_001 = 2.0 * range;
        Implicit _lattice_002 = Divide(_lattice_000, _lattice_001);
        Implicit lattice = Subtract(_lattice_002, bias);
        Implicit inverse = Multiply(-1.0, lattice);
        return inverse;
    }

    static Implicit thinLattice(Vector3d p) {
        double _tpms_000 = 1.0 - gyroid;
        Implicit x = new Implicit(p.x, new Vector3d(1.0, 0.0, 0.0));
        double halfsizeX = size_x * 0.5;
        double frequencyX = PI / halfsizeX;
        Implicit _sx_000 = Multiply(x, frequencyX);
        Implicit sx = Sin(_sx_000);
        Implicit _tpms_001 = Multiply(_tpms_000, sx);
        Implicit y = new Implicit(p.y, new Vector3d(0.0, 1.0, 0.0));
        double halfsizeY = size_y * 0.5;
        double frequencyY = PI / halfsizeY;
        Implicit _sy_000 = Multiply(y, frequencyY);
        Implicit sy = Sin(_sy_000);
        Implicit _tpms_002 = Multiply(_tpms_001, sy);
        Implicit z = new Implicit(p.z, new Vector3d(0.0, 0.0, 1.0));
        double halfsizeZ = size_z * 0.5;
        double frequencyZ = PI / halfsizeZ;
        Implicit _sz_000 = Multiply(z, frequencyZ);
        Implicit sz = Sin(_sz_000);
        Implicit _tpms_003 = Multiply(_tpms_002, sz);
        Implicit _tpms_004 = Multiply(drop_x, sx);
        Implicit _cy_000 = Multiply(y, frequencyY);
        Implicit cy = Cos(_cy_000);
        Implicit _tpms_005 = Multiply(_tpms_004, cy);
        Implicit _cz_000 = Multiply(z, frequencyZ);
        Implicit cz = Cos(_cz_000);
        Implicit _tpms_006 = mix(cz, 1.0, gyroid);
        Implicit _tpms_007 = Multiply(_tpms_005, _tpms_006);
        Implicit _tpms_008 = Add(_tpms_003, _tpms_007);
        Implicit _tpms_009 = Multiply(drop_y, sy);
        Implicit _tpms_010 = Multiply(_tpms_009, cz);
        Implicit _cx_000 = Multiply(x, frequencyX);
        Implicit cx = Cos(_cx_000);
        Implicit _tpms_011 = mix(cx, 1.0, gyroid);
        Implicit _tpms_012 = Multiply(_tpms_010, _tpms_011);
        Implicit _tpms_013 = Add(_tpms_008, _tpms_012);
        Implicit _tpms_014 = Multiply(drop_z, sz);
        Implicit _tpms_015 = Multiply(_tpms_014, cx);
        Implicit _tpms_016 = mix(cy, 1.0, gyroid);
        Implicit _tpms_017 = Multiply(_tpms_015, _tpms_016);
        Implicit tpms = Add(_tpms_013, _tpms_017);
        Implicit _lattice_000 = Multiply(tpms, size_x);
        double infimum = -1.0 * SQRT2;
        double range = supremum - infimum;
        double _lattice_001 = 2.0 * range;
        Implicit _lattice_002 = Divide(_lattice_000, _lattice_001);
        Implicit lattice = Subtract(_lattice_002, bias);
        Implicit _thin_000 = Abs(lattice);
        double _thin_001 = thickness * 0.5;
        Implicit thin = Subtract(_thin_000, _thin_001);
        return thin;
    }

    static Implicit twinLattice(Vector3d p) {
        Implicit thin = thinLattice(p);
        Implicit twin = Multiply(-1.0, thin);
        return twin;
    }

    public static Implicit IndexedLattice(Vector3d p) {
        switch (LatticeIndex) {
            case 0: return solidLattice(p);
            case 1: return inverseLattice(p);
            case 2: return thinLattice(p);
            case 3: return twinLattice(p);
        }
        return Sphere(p, new Vector3d(0.0), 5.0);
    }

    public static Implicit ScaledLattice(Vector3d scaledP) {
        Vector3d p = (scaledP - center) * 10.0;
        Implicit lattice = IndexedLattice(p);
        return Divide(lattice, 10.0);
    }

}
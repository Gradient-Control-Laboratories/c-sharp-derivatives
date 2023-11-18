using g4;

public class ImplicitFromCode : BoundedImplicitFunction3d{
    public string Code;
    public AxisAlignedBox3d BBox;

    public ImplicitFromCode(AxisAlignedBox3d bbox) {
        BBox = bbox;
    }

    public double Value(ref Vector3d p) {
        return IndexedLattice(p).Distance;
    }

    public AxisAlignedBox3d Bounds()
    {
        return BBox;
    }

    /////////////
    ///
    public static int VariantIndex;
    public static double size_x;
    public static double size_y;
    public static double size_z;
    public static double bias;
    public static double thickness;
    public static double length;
    public static double rotateA;
    public static double rotateB;

    public static double rotateADeriv;
    public static double rotateBDeriv;

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

    static Implicit latticeLattice(Vector3d p)
    {
        Implicit x = new Implicit(p.x, new Vector3d(1.0, 0.0, 0.0));
        Implicit y = new Implicit(p.y, new Vector3d(0.0, 1.0, 0.0));
        Implicit z = new Implicit(p.z, new Vector3d(0.0, 0.0, 1.0));
        double halfLength = length * 0.5;
        double _beam_000 = -1.0 * halfLength;
        Implicit beam = BoxCenter(p, Vector3d.Zero, new Vector3d(2.0, 2.0, 5.0));
        //Implicit beam = LineSegment(x, y, z, 0.0, 0.0, _beam_000, 0.0, 0.0, halfLength);
        Implicit lattice = Subtract(beam, bias);
        return lattice;
    }

    static Implicit rotatedLattice(Vector3d p)
    {
        // var rotation = Matrix3d.AxisAngleD(new Vector3d(1.0, 0.0, 0.0), rotateA * 180);  // API is in degrees
        var sA = Math.Sin(PI * rotateA);
        var cA = Math.Cos(PI * rotateA);
        var rotation = new Matrix3d(1, 0, 0, 0, cA, -sA, 0, sA, cA);
        var rotationDeriv = new Matrix3d(1, 0, 0, 0, -sA, -cA, 0, cA, -sA);

        var rotatedP = rotation * p;
        var rotatedPDeriv = rotationDeriv * p;
        var rotatedImplicit = latticeLattice(rotatedP);

        rotateADeriv = Vector3d.Dot(rotatedImplicit.Gradient, rotatedPDeriv);
        return new Implicit(rotateADeriv, rotatedImplicit.Gradient);
    //    rotateADeriv = latticeLattice(p).Gradient.Length;
        return rotatedImplicit;
    }

    public static Implicit IndexedLattice(Vector3d p)
    {
        Implicit lattice = rotatedLattice(p);
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


/// <summary>
/// 
/// </summary>


    const double PI = 3.14159265358979;
    const double SQRT2 = 1.41421356237;
    const double SQRT3 = 1.73205080757;

    static Implicit CreateImplicit() { return new Implicit(0.0, new Vector3d(0.0)); }
    static Implicit CreateImplicit(double iValue) { return new Implicit(iValue, new Vector3d(0.0)); }

    static double mix(double a, double b, double t) {
        return (1 - t) * a + t * b;
    }

    static Vector3d mix(Vector3d a, Vector3d b, double t) {
        return (1 - t) * a + t * b;
    }

    static Implicit mix(Implicit a, Implicit b, double t) {
        return new Implicit(
            mix(a.Distance, b.Distance, t),
            mix(a.Gradient, b.Gradient, t)
        );
    }

    static Implicit Negate(Implicit iImplicit) {
        return new Implicit(-iImplicit.Distance, -iImplicit.Gradient);
    }

    static Implicit Add(Implicit a, Implicit b) {
        return new Implicit(a.Distance + b.Distance, a.Gradient + b.Gradient);
    }

    static Implicit Subtract(Implicit a, Implicit b)  {
        return new Implicit(a.Distance - b.Distance, a.Gradient - b.Gradient);
    }

    static Implicit Add(double iT, Implicit iImplicit) {
        return new Implicit(iT + iImplicit.Distance, iImplicit.Gradient);
    }
    static Implicit Add(Implicit iImplicit, double iT) { return Add(iT, iImplicit); }
    static Implicit Subtract(double iT, Implicit iImplicit) { return Add(iT, Negate(iImplicit)); }
    static Implicit Subtract(Implicit iImplicit, double iT) { return Add(-iT, iImplicit); }

    static Implicit Multiply(Implicit a, Implicit b) {
        return new Implicit(a.Distance * b.Distance, a.Distance * b.Gradient + b.Distance * a.Gradient);
    }
    static Implicit Multiply(double iT, Implicit iImplicit) { return new Implicit(iT * iImplicit.Distance, iT * iImplicit.Gradient); }
    static Implicit Multiply(Implicit iImplicit, double iT) { return Multiply(iT, iImplicit); }

    static Implicit Divide(Implicit a, Implicit b) {
        return new Implicit(a.Distance / b.Distance, (b.Distance * a.Gradient - a.Distance * b.Gradient) / (b.Distance * b.Distance));
    }
    static Implicit Divide(Implicit a, double b) { return new Implicit(a.Distance / b, a.Gradient / b); }

    static Implicit Min(Implicit a, Implicit b) {
        if (a.Distance <= b.Distance)
            return a;

        return b;
    }
    static Implicit Min(Implicit a, Implicit b, Implicit c) { return Min(a, Min(b, c)); }
    static Implicit Min(Implicit a, double b) { return Min(a, CreateImplicit(b)); }
    static Implicit Min(double a, Implicit b) { return Min(b, CreateImplicit(a)); }
    static Implicit Min(Implicit a, Implicit b, Implicit c, Implicit d) { return Min(a, Min(b, Min(c, d))); }

    static Implicit Max(Implicit a, Implicit b) {
        if (a.Distance >= b.Distance)
            return a;

        return b;
    }
    static Implicit Max(Implicit a, Implicit b, Implicit c) { return Max(a, Max(b, c)); }
    static Implicit Max(Implicit a, double b) { return Max(a, CreateImplicit(b)); }
    static Implicit Max(double a, Implicit b) { return Max(b, CreateImplicit(a)); }
    static Implicit Max(Implicit a, Implicit b, Implicit c, Implicit d) { return Max(a, Max(b, Max(c, d))); }

    static Vector3d Max(Vector3d a, Vector3d b) {
        return new Vector3d(Math.Max(a.x, b.x), Math.Max(a.y, b.y), Math.Max(a.z, b.z));
    }

    static Implicit Compare(Implicit iA, Implicit iB) {
        if (iA.Distance < iB.Distance)
            return CreateImplicit(-1.0);
            
        if (iA.Distance > iB.Distance)
            return CreateImplicit(1.0);

        return CreateImplicit(0.0);
    }
    static Implicit Compare(Implicit iA, double iB) { return Compare(iA, CreateImplicit(iB)); }

    static Implicit Exp(Implicit iImplicit) {
        double exp = Math.Exp(iImplicit.Distance);
        return new Implicit(exp, exp * iImplicit.Gradient);
    }

    static Implicit Log(Implicit iImplicit) {
        return new Implicit(Math.Log(iImplicit.Distance), iImplicit.Gradient / iImplicit.Distance);
    }

    static Implicit Pow(Implicit iMantissa, Implicit iExponent) {
        double result = Math.Pow(iMantissa.Distance, iExponent.Distance);
        return new Implicit(result, result * Math.Log(iMantissa.Distance) * iMantissa.Gradient);
    }

    static Implicit Sqrt(Implicit iImplicit) {
        double sqrt = Math.Sqrt(iImplicit.Distance);
        return new Implicit(sqrt, iImplicit.Gradient / (2.0 * sqrt));
    }

    static Implicit Square(Implicit iImplicit) {
        double square = iImplicit.Distance * iImplicit.Distance ;
        return new Implicit(square, 2 * square * iImplicit.Gradient);
    }

    static Implicit Abs(Implicit iImplicit) {
        return new Implicit(Math.Abs(iImplicit.Distance), Math.Sign(iImplicit.Distance) * iImplicit.Gradient);
    }

    static Vector3d Abs(Vector3d v) {
        return new Vector3d(Math.Abs(v.x), Math.Abs(v.y), Math.Abs(v.z));
    }

    static Implicit Sin(Implicit iImplicit) {
        return new Implicit(Math.Sin(iImplicit.Distance), Math.Cos(iImplicit.Distance) * iImplicit.Gradient);
    }

    static Implicit Cos(Implicit iImplicit) {
        return new Implicit(Math.Cos(iImplicit.Distance), -Math.Sin(iImplicit.Distance) * iImplicit.Gradient);
    }

    static Implicit Tan(Implicit iImplicit) {
        double sec = 1 / Math.Cos(iImplicit.Distance);
        return new Implicit(Math.Tan(iImplicit.Distance), sec * sec * iImplicit.Gradient);
    }

    static Vector3d Boundary(Vector3d iP, Implicit i) {
        return -i.Distance * i.Gradient;
    }

    static Implicit Shell(Implicit iImplicit, double thickness, double bias) {
        thickness *= 0.5;
        return Subtract(Abs(Add(iImplicit, bias * thickness)), thickness);
    }

    static Implicit EuclideanNorm(Implicit a, Implicit b) {
        return Sqrt(Add(Multiply(a, a), Multiply(b, b)));
    }
    static Implicit EuclideanNorm(Implicit a, Implicit b, Implicit c) {
        return Sqrt(Add(Add(Multiply(a, a), Multiply(b, b)), Multiply(c, c)));
    }

    // Booleans

    // https://mercury.sexy/hg_sdf/
    static Implicit UnionEuclidean(Implicit a, Implicit b, double radius) {
        Implicit ua = Subtract(a, radius);
        Implicit ub = Subtract(b, radius);

        Implicit ab = Min(ua, ub);

        if (Math.Max(ua.Distance, ub.Distance) < 0) {
            ab = Negate(EuclideanNorm(ua, ub));
        }

        return Add(ab, radius);
    }

    static Implicit IntersectionEuclidean(Implicit a, Implicit b, double radius) {
        return Negate(UnionEuclidean(Negate(a), Negate(b), radius));
    }

    static Implicit IntersectionSmoothExp(Implicit a, Implicit b, double radius) {
        return Negate(UnionEuclidean(Negate(a), Negate(b), radius));
    }


    static Implicit UnionEuclidean(Implicit a, Implicit b, Implicit c, double radius) {
        Implicit zero = CreateImplicit(0.0);
        Implicit r = CreateImplicit(radius);
        Implicit ua = Max(Subtract(r, a), zero);
        Implicit ub = Max(Subtract(r, b), zero);
        Implicit uc = Max(Subtract(r, c), zero);

        Implicit abc = Min(a, Min(b, c));

        Implicit op = Subtract(Max(r, abc), EuclideanNorm(ua, ub, uc));

        if (abc.Distance > 0.0)
            op.Gradient = abc.Gradient;

        Vector3d bary = 0.33 * new Vector3d(ua.Distance, ub.Distance, uc.Distance) / (ua.Distance + ub.Distance + uc.Distance);

        return op;
    }

    static Implicit UnionChamfer(Implicit iA, Implicit iB, double k) {
        Implicit h = Multiply(Max(Subtract(CreateImplicit(k), Abs(Subtract(iA, iB))), CreateImplicit()), 1.0 / k);
        Implicit h2 = Multiply(h, 0.5);
        Implicit result = Subtract(Min(iA, iB), Multiply(h2, k * 0.5));
        double param = h2.Distance;

        return result;
    }

    static Implicit UnionRound(Implicit iA, Implicit iB, double k) {
        Implicit h = Multiply(Max(Subtract(CreateImplicit(k), Abs(Subtract(iA, iB))), CreateImplicit()), 1.0 / k);
        Implicit h2 = Multiply(Multiply(h, h), 0.5);
        Implicit result = Subtract(Min(iA, iB), Multiply(h2, k * 0.5));
        double param = h2.Distance;

        return result;
    }

    // Primitives

    static Implicit Plane(Vector3d p, Vector3d origin, Vector3d normal) {
        Vector3d grad = normal.Normalized;
        double v = Vector3d.Dot(p - origin, grad);
        return new Implicit(v, grad);
    }
    static Implicit Plane(Vector2d p, Vector2d origin, Vector2d normal) {
        return Plane(new Vector3d(p.x, p.y, 0.0), new Vector3d(origin.x, origin.y, 0.0), new Vector3d(normal.x, normal.y, 0.0));
    }


    static Implicit Circle(Vector2d p, Vector2d center, double iRadius) {
        Vector2d centered = p - center;
        double len = centered.Length;
        double length = len - iRadius;
        Vector2d normalized = centered / len;
        return new Implicit(length, new Vector3d(normalized.x, normalized.y, 0.0));
    }

    static Matrix2d Rotate2(double theta) {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        return new Matrix2d(
            new Vector2d(c, -s),
            new Vector2d(s, c)
        );
    }

    static Matrix3d RotateZ(double theta) {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        return new Matrix3d(
            new Vector3d(c, -s, 0.0),
            new Vector3d(s, c, 0.0),
            new Vector3d(0.0, 0.0, 1.0),
            true
        );
    }

    static Vector3d RotateX(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(p.x, -sa * p.y + ca * p.z, ca * p.y + sa * p.z);
    }
    static Vector3d RotateY(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(ca * p.x + sa * p.z, p.y, -sa * p.x + ca * p.z);
    }
    static Vector3d RotateZ(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(ca * p.x + sa * p.y, -sa * p.x + ca * p.y, p.z);
    }

    static Implicit Sphere(Vector3d iP, Vector3d iCenter, double iRadius) {
        Vector3d centered = iP - iCenter;
        double length = centered.Length;
        double dist = length - iRadius;
        return new Implicit(dist, centered / length);
    }

    static Implicit BoxCenter(Vector3d iP, Vector3d iCenter, Vector3d iSize) {
        Vector3d p = iP - iCenter;
        Vector3d b = iSize * 0.5;

        Vector3d d = Abs(p) - b;
        Vector3d maxvec = Max(d, Vector3d.Zero);
        double dist = maxvec.Length + Math.Min(Math.Max(d.x, Math.Max(d.y, d.z)), 0);

        Vector3d grad = (d.x > d.y) && (d.x > d.z) ? new Vector3d(1, 0, 0) :
            (d.y > d.z ? new Vector3d(0, 1, 0) : new Vector3d(0, 0, 1));

        if (d.x > 0 || d.y > 0 || d.z > 0)
        {
            d = Max(d, Vector3d.Zero);
            grad = d / d.Length;
        }

        grad *= new Vector3d(Math.Sign(p.x), Math.Sign(p.y), Math.Sign(p.z));
        return new Implicit(dist, grad);
    }

    static Implicit Box(Vector3d iP, Vector3d iA, Vector3d iB) {
        Vector3d center = (iA + iB) * 0.5;
        Vector3d size = Abs(iA - iB);
        return BoxCenter(iP, center, size);
    }

    static Vector3d center = new Vector3d(0.5);

}
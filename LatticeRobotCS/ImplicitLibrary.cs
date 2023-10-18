/////////////////////////////////////////////////////////////////////////////////////////////
//
// LatticeRobot CodeRep output v0.0
// This code in the 'Common' tab is test output from LatticeRobot's lattice CodeRep generator.  
// It is intended for other implicits to use to fill shapes with spatially varying lattices.
//
// https://www.latticerobot.com
//
// LICENSE: This code is provided under an MIT license.  
// BETA WARNING: This format will evolve through the fall of 2023
//
// Documentation:
// Immediately below, there is a library, followed by the custom code for this parameterized
// unit cell (line 336).  A complier directive replaces the uniforms with internal parameters
// for the purpose of this Shadertoy demo.
//
// Lattice Variants and Index:
// The latticeIndex uniform accesses the four variants of the lattice F(p), with the shape 
// defined by the preimage of the non-positive values of that variant:
//
//   0: Solid. The Solid lattice:
//       F(p) - bias
//
//   1: Inverse. The inverse of the Solid lattice:
//       -F(p) + bias
//
//   2: Thin. A thin band of the Solid lattice:
//       abs(F(p) - bias) - thickness / 2
//
//   3: Twin. The twin axes produced by the inverse of the Thin lattice:
//       -abs(F(p) - bias) + thickness / 2
//
// Common Lattice Parameters (with default value and units):
//   * size_x, size_y, and size_z (10 mm): The dimensions of the desired bounding box. 
//   * bias (0 mm): The offset of the baseline lattice surface.  Not precise distance.
//   * thickness (1 mm): Thin and Twin only.  The thickness of the wall of the thin lattice or
//        the span of the spacing of the twin lattice.  Not precise distance.
//
// Special Lattice Parameters (for only this lattice):
//   * drop_x, drop_y, and drop_z (1): Lattice-specific parameters that attenuate coefficients.
//   * gyroid (0): interpolates the Diamond TPMS to Gyroid TPMS.  
//  
/////////////////////////////////////////////////////////////////////////////////////////////

using g4;
public partial class DiamondImplicit {

    const double PI = 3.14159265358979;
    const double SQRT2 = 1.41421356237;
    const double SQRT3 = 1.73205080757;

    Implicit CreateImplicit() { return new Implicit(0.0, new Vector3d(0.0)); }
    Implicit CreateImplicit(double iValue) { return new Implicit(iValue, new Vector3d(0.0)); }

    double mix(double a, double b, double t) {
        return (1 - t) * a + t * b;
    }

    Vector3d mix(Vector3d a, Vector3d b, double t) {
        return (1 - t) * a + t * b;
    }

    Implicit mix(Implicit a, Implicit b, double t) {
        return new Implicit(
            mix(a.Distance, b.Distance, t),
            mix(a.Gradient, b.Gradient, t)
        );
    }

    Implicit Negate(Implicit iImplicit) {
        return new Implicit(-iImplicit.Distance, -iImplicit.Gradient);
    }

    Implicit Add(Implicit a, Implicit b) {
        return new Implicit(a.Distance + b.Distance, a.Gradient + b.Gradient);
    }

    Implicit Subtract(Implicit a, Implicit b)  {
        return new Implicit(a.Distance - b.Distance, a.Gradient - b.Gradient);
    }

    Implicit Add(double iT, Implicit iImplicit) {
        return new Implicit(iT + iImplicit.Distance, iImplicit.Gradient);
    }
    Implicit Add(Implicit iImplicit, double iT) { return Add(iT, iImplicit); }
    Implicit Subtract(double iT, Implicit iImplicit) { return Add(iT, Negate(iImplicit)); }
    Implicit Subtract(Implicit iImplicit, double iT) { return Add(-iT, iImplicit); }

    Implicit Multiply(Implicit a, Implicit b) {
        return new Implicit(a.Distance * b.Distance, a.Distance * b.Gradient + b.Distance * a.Gradient);
    }
    Implicit Multiply(double iT, Implicit iImplicit) { return new Implicit(iT * iImplicit.Distance, iT * iImplicit.Gradient); }
    Implicit Multiply(Implicit iImplicit, double iT) { return Multiply(iT, iImplicit); }

    Implicit Divide(Implicit a, Implicit b) {
        return new Implicit(a.Distance / b.Distance, (b.Distance * a.Gradient - a.Distance * b.Gradient) / (b.Distance * b.Distance));
    }
    Implicit Divide(Implicit a, double b) { return new Implicit(a.Distance / b, a.Gradient / b); }

    Implicit Min(Implicit a, Implicit b) {
        if (a.Distance <= b.Distance)
            return a;

        return b;
    }
    Implicit Min(Implicit a, Implicit b, Implicit c) { return Min(a, Min(b, c)); }
    Implicit Min(Implicit a, Implicit b, Implicit c, Implicit d) { return Min(a, Min(b, Min(c, d))); }

    Implicit Max(Implicit a, Implicit b) {
        if (a.Distance >= b.Distance)
            return a;

        return b;
    }
    Implicit Max(Implicit a, Implicit b, Implicit c) { return Max(a, Max(b, c)); }
    Implicit Max(Implicit a, double b) { return Max(a, CreateImplicit(b)); }
    Implicit Max(Implicit a, Implicit b, Implicit c, Implicit d) { return Max(a, Max(b, Max(c, d))); }

    Vector3d Max(Vector3d a, Vector3d b) {
        return new Vector3d(Math.Max(a.x, b.x), Math.Max(a.y, b.y), Math.Max(a.z, b.z));
    }

    Implicit Compare(Implicit iA, Implicit iB) {
        if (iA.Distance < iB.Distance)
            return CreateImplicit(-1.0);
            
        if (iA.Distance > iB.Distance)
            return CreateImplicit(1.0);

        return CreateImplicit(0.0);
    }
    Implicit Compare(Implicit iA, double iB) { return Compare(iA, CreateImplicit(iB)); }

    Implicit Exp(Implicit iImplicit) {
        double exp = Math.Exp(iImplicit.Distance);
        return new Implicit(exp, exp * iImplicit.Gradient);
    }

    Implicit Log(Implicit iImplicit) {
        return new Implicit(Math.Log(iImplicit.Distance), iImplicit.Gradient / iImplicit.Distance);
    }

    Implicit Pow(Implicit iMantissa, Implicit iExponent) {
        double result = Math.Pow(iMantissa.Distance, iExponent.Distance);
        return new Implicit(result, result * Math.Log(iMantissa.Distance) * iMantissa.Gradient);
    }

    Implicit Sqrt(Implicit iImplicit) {
        double sqrt = Math.Sqrt(iImplicit.Distance);
        return new Implicit(sqrt, iImplicit.Gradient / (2.0 * sqrt));
    }

    Implicit Square(Implicit iImplicit) {
        double square = iImplicit.Distance * iImplicit.Distance ;
        return new Implicit(square, 2 * square * iImplicit.Gradient);
    }

    Implicit Abs(Implicit iImplicit) {
        return new Implicit(Math.Abs(iImplicit.Distance), Math.Sign(iImplicit.Distance) * iImplicit.Gradient);
    }

    Vector3d Abs(Vector3d v) {
        return new Vector3d(Math.Abs(v.x), Math.Abs(v.y), Math.Abs(v.z));
    }

    Implicit Sin(Implicit iImplicit) {
        return new Implicit(Math.Sin(iImplicit.Distance), Math.Cos(iImplicit.Distance) * iImplicit.Gradient);
    }

    Implicit Cos(Implicit iImplicit) {
        return new Implicit(Math.Cos(iImplicit.Distance), -Math.Sin(iImplicit.Distance) * iImplicit.Gradient);
    }

    Implicit Tan(Implicit iImplicit) {
        double sec = 1 / Math.Cos(iImplicit.Distance);
        return new Implicit(Math.Tan(iImplicit.Distance), sec * sec * iImplicit.Gradient);
    }

    Vector3d Boundary(Vector3d iP, Implicit i) {
        return -i.Distance * i.Gradient;
    }

    Implicit Shell(Implicit iImplicit, double thickness, double bias) {
        thickness *= 0.5;
        return Subtract(Abs(Add(iImplicit, bias * thickness)), thickness);
    }

    Implicit EuclideanNorm(Implicit a, Implicit b) {
        return Sqrt(Add(Multiply(a, a), Multiply(b, b)));
    }
    Implicit EuclideanNorm(Implicit a, Implicit b, Implicit c) {
        return Sqrt(Add(Add(Multiply(a, a), Multiply(b, b)), Multiply(c, c)));
    }

    // Booleans

    // https://mercury.sexy/hg_sdf/
    Implicit UnionEuclidean(Implicit a, Implicit b, double radius) {
        Implicit ua = Subtract(a, radius);
        Implicit ub = Subtract(b, radius);

        Implicit ab = Min(ua, ub);

        if (Math.Max(ua.Distance, ub.Distance) < 0) {
            ab = Negate(EuclideanNorm(ua, ub));
        }

        return Add(ab, radius);
    }

    Implicit IntersectionEuclidean(Implicit a, Implicit b, double radius) {
        return Negate(UnionEuclidean(Negate(a), Negate(b), radius));
    }

    Implicit IntersectionSmoothExp(Implicit a, Implicit b, double radius) {
        return Negate(UnionEuclidean(Negate(a), Negate(b), radius));
    }


    Implicit UnionEuclidean(Implicit a, Implicit b, Implicit c, double radius) {
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

    Implicit UnionChamfer(Implicit iA, Implicit iB, double k) {
        Implicit h = Multiply(Max(Subtract(CreateImplicit(k), Abs(Subtract(iA, iB))), CreateImplicit()), 1.0 / k);
        Implicit h2 = Multiply(h, 0.5);
        Implicit result = Subtract(Min(iA, iB), Multiply(h2, k * 0.5));
        double param = h2.Distance;

        return result;
    }

    Implicit UnionRound(Implicit iA, Implicit iB, double k) {
        Implicit h = Multiply(Max(Subtract(CreateImplicit(k), Abs(Subtract(iA, iB))), CreateImplicit()), 1.0 / k);
        Implicit h2 = Multiply(Multiply(h, h), 0.5);
        Implicit result = Subtract(Min(iA, iB), Multiply(h2, k * 0.5));
        double param = h2.Distance;

        return result;
    }

    // Primitives

    Implicit Plane(Vector3d p, Vector3d origin, Vector3d normal) {
        Vector3d grad = normal.Normalized;
        double v = Vector3d.Dot(p - origin, grad);
        return new Implicit(v, grad);
    }
    Implicit Plane(Vector2d p, Vector2d origin, Vector2d normal) {
        return Plane(new Vector3d(p.x, p.y, 0.0), new Vector3d(origin.x, origin.y, 0.0), new Vector3d(normal.x, normal.y, 0.0));
    }


    Implicit Circle(Vector2d p, Vector2d center, double iRadius) {
        Vector2d centered = p - center;
        double len = centered.Length;
        double length = len - iRadius;
        Vector2d normalized = centered / len;
        return new Implicit(length, new Vector3d(normalized.x, normalized.y, 0.0));
    }

    Matrix2d Rotate2(double theta) {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        return new Matrix2d(
            new Vector2d(c, -s),
            new Vector2d(s, c)
        );
    }

    Matrix3d RotateZ(double theta) {
        double c = Math.Cos(theta);
        double s = Math.Sin(theta);
        return new Matrix3d(
            new Vector3d(c, -s, 0.0),
            new Vector3d(s, c, 0.0),
            new Vector3d(0.0, 0.0, 1.0),
            true
        );
    }

    Vector3d RotateX(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(p.x, -sa * p.y + ca * p.z, ca * p.y + sa * p.z);
    }
    Vector3d RotateY(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(ca * p.x + sa * p.z, p.y, -sa * p.x + ca * p.z);
    }
    Vector3d RotateZ(Vector3d p, double a) {
        double sa = Math.Sin(a);
        double ca = Math.Cos(a);
        return new Vector3d(ca * p.x + sa * p.y, -sa * p.x + ca * p.y, p.z);
    }

    Implicit Sphere(Vector3d iP, Vector3d iCenter, double iRadius) {
        Vector3d centered = iP - iCenter;
        double length = centered.Length;
        double dist = length - iRadius;
        return new Implicit(dist, centered / length);
    }

    Implicit BoxCenter(Vector3d iP, Vector3d iCenter, Vector3d iSize) {
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

    Implicit Box(Vector3d iP, Vector3d iA, Vector3d iB) {
        Vector3d center = (iA + iB) * 0.5;
        Vector3d size = Abs(iA - iB);
        return BoxCenter(iP, center, size);
    }

    Vector3d center = new Vector3d(0.5);

}
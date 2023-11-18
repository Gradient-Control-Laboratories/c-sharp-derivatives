using g4;
using System;

public struct Implicit {
    double distance = 0;
    public double Distance {
        get { return distance; }
        set { distance = value;
            //distance = Double.IsNaN(value) ? 0 : value; 
            // if (Double.IsNaN(value))
            //     distance = value;
            // else
            //     distance = 0;
            
            }
    }
    public Vector3d Gradient { get; set; }

    public Implicit(double distance, Vector3d gradient) {
        this.distance = Double.IsNaN(distance) ? 0 : distance;;
        Gradient = gradient;
    }
}
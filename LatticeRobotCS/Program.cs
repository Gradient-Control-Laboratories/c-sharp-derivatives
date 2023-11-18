﻿using System;
// using System.CodeDom;
// using System.CodeDom.Compiler;

using g4;
using Microsoft.VisualBasic;

internal class Program {
    private static int Main(string[] args) {

        // Not in use - to dynamically load a CodeRep
        // var unitCellLocation = "Test-Rotating-Shape";

        // if (args.Length == 0) 
        //     Console.WriteLine($"Using default unit cell: {unitCellLocation}.");
        // else
        //     unitCellLocation = args[0];
        
        // var unitCell = new ImplicitUnitCell(Path.Combine(@"..\..\..\..\Samples\", unitCellLocation), 2);
        // With this implementation of ImplicitUnitCell, we can only set constant parameters.  
        // unitCell.SetParameter("thickness", 2);
        // unitCell.SetParameter("length", 5);
        // unitCell.VariantIndex = (int)LatticeVariant.thin;

        Console.WriteLine("Computing mesh...");

        // generateMeshF() meshes the input implicit function at
        // the given cell resolution, and writes out the resulting mesh    
        Action<BoundedImplicitFunction3d, int, string> generateMeshF = (root, numcells, path) => {
            MarchingCubes c = new MarchingCubes();
            c.Implicit = root;
            c.RootMode = MarchingCubes.RootfindingModes.LerpSteps;      // cube-edge convergence method
            c.RootModeSteps = 5;                                        // number of iterations
            c.Bounds = root.Bounds();
            c.CubeSize = c.Bounds.MaxDim / numcells;
            c.Bounds.Expand(3 * c.CubeSize);                            // leave a buffer of cells
            c.Generate();
            MeshNormals.QuickCompute(c.Mesh);                           // generate normals
            StandardMeshWriter.WriteMesh(path, c.Mesh, WriteOptions.Defaults);   // write mesh
        };

        AxisAlignedBox3d bbox = new(-5 * Vector3d.One, 5 * Vector3d.One);
        ImplicitBox3d box = new ImplicitBox3d() { Box = new Box3d(bbox) };
        var unitCell = new ImplicitFromCode(bbox);

        generateMeshF(new ImplicitIntersection3d() { A = unitCell, B = box }, 128, "LatticeRobot.obj");

        return 0;
    }

}
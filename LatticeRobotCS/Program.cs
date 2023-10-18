using g4;

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

ImplicitBox3d box = new ImplicitBox3d() {
    Box = new Box3d(bbox)
};

ImplicitUnitCell unitcell = new(2);

Console.WriteLine("Computing mesh...");
generateMeshF(new ImplicitIntersection3d() { A = unitcell, B = box }, 128, "LatticeRobot.obj");


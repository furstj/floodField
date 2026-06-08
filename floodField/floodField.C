/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    floodField

Description
    Sets initial alpha.water field for interFoam simulations with free surface
    flow. Uses a flood-fill algorithm starting from seed points to fill all
    cells below the water level, excluding solid regions (e.g., empty boat hull).

    The water surface is defined by a point on the surface and the gravity
    direction, allowing arbitrary coordinate system orientations.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "DynamicList.H"
#include "IFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Sets initial alpha.water field for interFoam using flood-fill\n"
        "from seed points, filling all cells below the water surface."
    );

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    const word dictName("floodFieldDict");

    IOdictionary floodFieldDict
    (
        IOobject
        (
            dictName,
            runTime.system(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Read gravity vector from constant/g
    fileName gFilePath = runTime.constant() / "g";
    IFstream gStream(gFilePath);
    dictionary gDict(gStream());

    vector gValue = gDict.get<vector>("value");
    dimensionedVector g("g", dimAcceleration, gValue);
    Info<< "Gravity vector: " << g.value() << endl;

    // Direction opposite to gravity (upward direction)
    vector negGDir = -g.value() / mag(g.value());
    Info<< "Upward direction (-g): " << negGDir << endl;

    // Read water level (scalar value representing projection along -g direction)
    scalar waterLevel = readScalar(floodFieldDict.lookup("waterLevel"));
    Info<< "Water level (projection value): " << waterLevel << endl;

    // Read seed points
    List<point> seedPoints(0);
    if (floodFieldDict.found("seedPoints"))
    {
        seedPoints = floodFieldDict.lookup("seedPoints");
    }

    Info<< "Number of seed points: " << seedPoints.size() << endl;

    // Create alpha.water field initialized to 0
    volScalarField alpha
    (
        IOobject
        (
            "alpha.water",
            runTime.timeName(),
            mesh,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedScalar("alpha0", dimless, 0.0)
    );

    Info<< "Initializing alpha.water field to 0..." << endl;
    scalarField& alphaF = alpha.ref();

    // Find cells containing seed points and add to queue
    DynamicList<label> queue(100);

    forAll(seedPoints, i)
    {
        const point& pt = seedPoints[i];
        label celli = mesh.findCell(pt);

        if (celli >= 0)
        {
            Info<< "Seed point " << pt << " -> cell " << celli << endl;
            alphaF[celli] = 1.0;
            queue.append(celli);
        }
        else
        {
            WarningInFunction
                << "Seed point " << pt << " not found in mesh, ignoring"
                << endl;
        }
    }

    if (queue.empty())
    {
        FatalErrorInFunction
            << "No valid seed points found. Cannot proceed with flood-fill."
            << exit(FatalError);
    }

    // Get cell centers and cell-cell connectivity
    const vectorField& C = mesh.C();
    const labelListList& cellCells = mesh.cellCells();

    // Track which cells are already in queue
    boolList inQueue(mesh.nCells(), false);
    forAll(queue, i)
    {
        inQueue[queue[i]] = true;
    }

    // Flood-fill algorithm (BFS)
    label filledCells = 0;

    while (!queue.empty())
    {
        label currentCell = queue.last();
        queue.pop_back();
        inQueue[currentCell] = false;
        filledCells++;

        // Check all neighbors
        forAll(cellCells[currentCell], j)
        {
            label neighborCell = cellCells[currentCell][j];

            // Only process if not already water and not in queue
            if (alphaF[neighborCell] == 0 && !inQueue[neighborCell])
            {
                // Check if neighbor is below water surface
                // Cell is below water if its projection along -g direction
                // is less than the water level
                scalar cellProjection = C[neighborCell] & negGDir;
                if (cellProjection < waterLevel)
                {
                    alphaF[neighborCell] = 1.0;
                    inQueue[neighborCell] = true;
                    queue.append(neighborCell);
                }
            }
        }
    }

    Info<< "Flood-fill complete. Filled " << filledCells << " cells with water."
        << endl;

    // Write the field
    alpha.write();

    Info<< "End" << endl;

    return 0;
}

// ************************************************************************* //

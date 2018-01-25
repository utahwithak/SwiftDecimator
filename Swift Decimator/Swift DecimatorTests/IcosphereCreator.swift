//
//  IcosphereCreator.swift
//  Leaf Press iOS
//
//  Created by Carl Wieland on 1/23/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation
import SceneKit
import simd
@testable import Swift_Decimator

public class IcoSphereCreator
{

    public static func icospherePoints(level: Int, radius: Float = 1, shift: float3 = float3(0,0,0)) -> ([float3], [Triangle]) {
        var positions = [float3]()
        var middlePointIndexCache = [Int: Int]()

        let addVertex: (_ vert: float3) -> Int = { vert in

            positions.append(((vert / sqrt(simd_dot(vert, vert))) * radius) + shift)
            return positions.count - 1
        }

        // return index of point in the middle of p1 and p2
        let getMiddlePoint: (Int, Int) -> Int = { p1, p2 in
            // first check if we have it already
            let firstIsSmaller = p1 < p2;
            let smallerIndex = firstIsSmaller ? p1 : p2;
            let greaterIndex = firstIsSmaller ? p2 : p1;
            let key = (smallerIndex << 32) + greaterIndex;


            if let val = middlePointIndexCache[key] {
                return val
            }

            // not in cache, calculate it
            let point1 = positions[p1];
            let point2 = positions[p2];
            let middle = (point1 + point2) / 2

            // add vertex makes sure point is on unit sphere
            let i = addVertex(middle);

            // store it, return index
            middlePointIndexCache[key] = i
            return i;
        }



        // create 12 vertices of a icosahedron
        let t = Float(1.0 + sqrt(5.0)) / 2.0;
        let initial = [float3(-1,  t,  0),
                     float3( 1,  t,  0),
                     float3(-1, -t,  0),
                     float3( 1, -t,  0),

                     float3( 0, -1,  t),
                     float3( 0,  1,  t),
                     float3( 0, -1, -t),
                     float3( 0,  1, -t),

                     float3( t,  0, -1),
                     float3( t,  0,  1),
                     float3(-t,  0, -1),
                     float3(-t,  0,  1)]
        initial.forEach({ _ = addVertex($0)})

        // create 20 triangles of the icosahedron
        var faces = [Triangle(v0: 0, v1: 11, v2: 5),
                     Triangle(v0: 0, v1: 5, v2: 1),
                     Triangle(v0: 0, v1: 1, v2: 7),
                     Triangle(v0: 0, v1: 7, v2: 10),
                     Triangle(v0: 0, v1: 10, v2: 11),

                     // 5 adjacent faces
            Triangle(v0: 1, v1: 5, v2: 9),
            Triangle(v0: 5, v1: 11, v2: 4),
            Triangle(v0: 11, v1: 10, v2: 2),
            Triangle(v0: 10, v1: 7, v2: 6),
            Triangle(v0: 7, v1: 1, v2: 8),

            // 5 faces around point 3
            Triangle(v0: 3, v1: 9, v2: 4),
            Triangle(v0: 3, v1: 4, v2: 2),
            Triangle(v0: 3, v1: 2, v2: 6),
            Triangle(v0: 3, v1: 6, v2: 8),
            Triangle(v0: 3, v1: 8, v2: 9),

            // 5 adjacent faces
            Triangle(v0: 4, v1: 9, v2: 5),
            Triangle(v0: 2, v1: 4, v2: 11),
            Triangle(v0: 6, v1: 2, v2: 10),
            Triangle(v0: 8, v1: 6, v2: 7),
            Triangle(v0: 9, v1: 8, v2: 1)]


        // refine triangles
        for _ in 0..<level {
            var faces2 = [Triangle]();
            for tri in faces {
                // replace triangle by 4 triangles
                let a = getMiddlePoint(tri.v0,tri.v1);
                let b = getMiddlePoint(tri.v1,tri.v2);
                let c = getMiddlePoint(tri.v2,tri.v0);

                faces2.append(Triangle(v0: tri.v0, v1: a, v2: c));
                faces2.append(Triangle(v0: tri.v1, v1: b, v2: a));
                faces2.append(Triangle(v0: tri.v2, v1: c, v2: b));
                faces2.append(Triangle(v0: a, v1: b, v2: c));
            }
            faces = faces2;
        }
        return (positions, faces)
    }

    static func normal(of a: float3, b: float3, c: float3) -> float3 {
        return simd_cross(b - a, c - a).normalized
    }
}

extension float3 {
    var normalized: float3 {
        return self / sqrt(simd_dot(self, self))
    }
}

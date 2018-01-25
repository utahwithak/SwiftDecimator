//
//  Decimator.swift
//  Swift Decimator
//
//  Created by Carl Wieland on 1/24/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation

struct Triangle {
    var v0 = 0
    var v1 = 0
    var v2 = 0
}

class Decimator {


    private var triangles: [Triangle]
    private var points: [Vector3]
    public var pointCount: Int {
        return points.count
    }

    public private(set) var initialTriangles = 0
    public private(set) var vertCount = 0
    public private(set) var triangleCount = 0
    public private(set) var edgeCount = 0
    private var diagBB = 0.0
    var vertices: [Vertex]
    var edges: [Edge]

    private var pqueue = PriorityQueue<EdgePriorityQueue>()

    private var trianglesTags = [Bool]()
    public var ecolManifoldConstraint = false


    init( points: [Vector3], triangles: [Triangle]) {
        self.points = points
        self.triangles = triangles
        initialTriangles = triangles.count
        triangleCount = triangles.count
        vertCount = points.count
        edges = [Edge](repeating:Edge(), count: triangleCount * 3)
        vertices = [Vertex](repeating:Vertex(), count: vertCount)

        var edge = Edge()
        edge.tag = true;
        edge.onBoundary = true;

        for t in 0..<triangleCount {
            let tri = [triangles[t].v0, triangles[t].v1, triangles[t].v2]
            trianglesTags.append(true)
            for k in 0..<3 {
                edge.v1 = tri[k];
                edge.v2 = tri[(k+1)%3];
                vertices[edge.v1].triangles.append(t)
                if let idEdge = getEdge(v1: edge.v1, v2: edge.v2) {
                    edges[idEdge].onBoundary = false;
                } else {
                    edges.append(edge);
                    vertices[edge.v1].edges.append(edgeCount);
                    vertices[edge.v2].edges.append(edgeCount);
                    edgeCount += 1
                }
            }
        }
        for v in 0..<vertCount {
            vertices[v].onBoundary = false;
            for itE in 0..<vertices[v].edges.count {
                let idEdge = vertices[v].edges[itE]
                if edges[idEdge].onBoundary {
                    vertices[v].onBoundary = true;
                    break;
                }
            }
        }


    }

    private func edgeCollapse(v1: Int, v2: Int) {

        for idTriangle in vertices[v2].triangles {
            var shift = 0
            var u = 0
            var w = 0
            if triangles[idTriangle].v0 == v2 {
                shift = 0;
                u = triangles[idTriangle].v1
                w = triangles[idTriangle].v2
            } else if triangles[idTriangle].v1 == v2 {
                shift = 1;
                u = triangles[idTriangle].v0
                w = triangles[idTriangle].v2
            } else {
                shift = 2;
                u = triangles[idTriangle].v0
                w = triangles[idTriangle].v1
            }

            if (u == v1) || (w == v1) {
                trianglesTags[idTriangle] = false;
                vertices[u].remove(triangle: idTriangle)
                vertices[w].remove(triangle: idTriangle)
                triangleCount -= 1
            }
            else if getTriangle(v1: v1, v2: u, v3: w) == -1 {
                vertices[v1].triangles.append(idTriangle)
                switch shift {
                case 0:
                    triangles[idTriangle].v0 = v1
                case 1:
                    triangles[idTriangle].v1 = v1
                default:
                    triangles[idTriangle].v2 = v1
                }
            } else {
                trianglesTags[idTriangle] = false;
                vertices[u].remove(triangle: idTriangle)
                vertices[w].remove(triangle: idTriangle)
                triangleCount -= 1
            }
        }
        var idEdge = 0
        for it in 0..<vertices[v2].edges.count  {
            if it >= vertices[v2].edges.count {
                break
            }
            idEdge = vertices[v2].edges[it]
            let w = (edges[idEdge].v1 == v2) ? edges[idEdge].v2 : edges[idEdge].v1;
            if w == v1 {
                edges[idEdge].tag = false;
                vertices[w].remove(edge: idEdge)
                edgeCount -= 1
            } else if getEdge(v1: v1, v2: w) == -1 {
                if edges[idEdge].v1 == v2 {
                    edges[idEdge].v1 = v1;
                } else{
                    edges[idEdge].v2 = v1;
                }
                vertices[v1].edges.append(idEdge);
            } else {
                edges[idEdge].tag = false;
                vertices[w].remove(edge: idEdge)
                edgeCount -= 1
            }
        }
        vertices[v2].tag = false;
        vertCount -= 1
        // update boundary edges
        var incidentVertices = [v1]
        for itE in 0..<vertices[v1].edges.count {
            incidentVertices.append((edges[idEdge].v1 != v1) ? edges[idEdge].v1 : edges[idEdge].v2)
            idEdge = vertices[v1].edges[itE]
            edges[idEdge].onBoundary = isBoundaryEdge(v1: edges[idEdge].v1, v2: edges[idEdge].v2) != -1
        }
        // update boundary vertices

        for idVertex in incidentVertices {
            vertices[idVertex].onBoundary = false;
            if vertices[idVertex].edges.contains(where: { edges[$0].onBoundary}) {
                vertices[idVertex].onBoundary = true
            }
        }
    }

    private func getTriangle(v1: Int, v2: Int, v3: Int) -> Int {

        for idTriangle in vertices[v1].triangles {
            let i = triangles[idTriangle].v0
            let j = triangles[idTriangle].v1
            let k = triangles[idTriangle].v2
            if ((i == v1 && j == v2 && k == v3) || (i == v1 && j == v3 && k == v2) ||
                (i == v2 && j == v1 && k == v3) || (i == v2 && j == v3 && k == v1) ||
                (i == v3 && j == v2 && k == v1) || (i == v3 && j == v1 && k == v2) )
            {
                return idTriangle;
            }
        }
        return -1;
    }

    private func getEdge(v1: Int,v2: Int) -> Int? {
        return vertices[v1].edges.first(where: { (idEdge) -> Bool in
            return (edges[idEdge].v1==v1 && edges[idEdge].v2==v2) || (edges[idEdge].v1==v2 && edges[idEdge].v2==v1)
        })
    }

    private func isBoundaryEdge(v1: Int, v2: Int) -> Int {
        var commonTri = -1;

        for itTriangle1 in vertices[v1].triangles {
            for itTriangle2 in vertices[v2].triangles {
                if (itTriangle1 == itTriangle2) {
                    if commonTri == -1 {
                        commonTri = itTriangle1;
                    } else {
                        return -1;
                    }
                }
            }
        }
        return commonTri;
    }

    private func isBoundaryVertex(v: Int) -> Bool {
        return vertices[v].edges.contains(where: { isBoundaryEdge(v1: edges[$0].v1, v2: edges[$0].v2) != -1 })
    }

    private func initializePriorityQueue() {
        var pqEdge = EdgePriorityQueue()

        for e in 0..<edges.count where edges[e].tag {
            let v1 = edges[e].v1;
            let v2 = edges[e].v2;
            if  (!ecolManifoldConstraint) || (manifoldConstraint(v1: v1, v2: v2)) {
                edges[e].qem = computeEdgeCost(v1: v1, v2: v2, newPos: &edges[e].position);
                pqEdge.qem = edges[e].qem
                pqEdge.name = e
                pqueue.push(pqEdge);
            }

        }
    }

    private func initializeQEM() {

        var coordMin = points[0];
        var coordMax = points[0];

        for p in 1..<points.count {
            let coord = points[p];
            if (coordMin.x > coord.x){ coordMin.x = coord.x;}
            if (coordMin.y > coord.y){ coordMin.y = coord.y;}
            if (coordMin.z > coord.z){ coordMin.z = coord.z;}
            if (coordMax.x < coord.x){ coordMax.x = coord.x;}
            if (coordMax.y < coord.y){ coordMax.y = coord.y;}
            if (coordMax.z < coord.z){ coordMax.z = coord.z;}
        }
        coordMax -= coordMin;
        diagBB = coordMax.magnitude

        var d = 0.0
        var area = 0.0
        for v in 0..<pointCount {
            vertices[v].resetQ()
            for idTriangle in vertices[v].triangles {
                let i = triangles[idTriangle].v0
                let j = triangles[idTriangle].v1
                let k = triangles[idTriangle].v2
                var n = (points[j] - points[i])^(points[k] - points[i]);
                area = n.magnitude
                n.normalize()
                d = -(points[v] * n)
                vertices[v].Q[0] += area * (n.x * n.x);
                vertices[v].Q[1] += area * (n.x * n.y);
                vertices[v].Q[2] += area * (n.x * n.z);
                vertices[v].Q[3] += area * (n.x * d);
                vertices[v].Q[4] += area * (n.y * n.y);
                vertices[v].Q[5] += area * (n.y * n.z);
                vertices[v].Q[6] += area * (n.y * d);
                vertices[v].Q[7] += area * (n.z * n.z);
                vertices[v].Q[8] += area * (n.z * d);
                vertices[v].Q[9] += area * (d     * d);
            }
        }

        let w = 1000.0
        for e in 0..<edges.count {
            let v1 = edges[e].v1
            let v2 = edges[e].v2
            let t = isBoundaryEdge(v1: v1, v2: v2)
            if t != -1 {
                let v3: Int
                if triangles[t].v0 != v1 && triangles[t].v0 != v2 {
                    v3 = triangles[t].v0
                } else if triangles[t].v1 != v1 && triangles[t].v1 != v2 {
                    v3 = triangles[t].v1
                } else {
                    v3 = triangles[t].v2
                }
                var u1 = points[v2] - points[v1];
                let u2 = points[v3] - points[v1];
                area = w * (u1^u2).magnitude
                u1.normalize()
                let n = (u2 - (u2 * u1) * u1).normalized


                d = -(points[v1] * n);
                vertices[v1].Q[0] += area * (n.x * n.x);
                vertices[v1].Q[1] += area * (n.x * n.y);
                vertices[v1].Q[2] += area * (n.x * n.z);
                vertices[v1].Q[3] += area * (n.x * d);
                vertices[v1].Q[4] += area * (n.y * n.y);
                vertices[v1].Q[5] += area * (n.y * n.z);
                vertices[v1].Q[6] += area * (n.y * d);
                vertices[v1].Q[7] += area * (n.z * n.z);
                vertices[v1].Q[8] += area * (n.z * d);
                vertices[v1].Q[9] += area * (d * d);

                d = -(points[v2] * n);
                vertices[v2].Q[0] += area * (n.x * n.x);
                vertices[v2].Q[1] += area * (n.x * n.y);
                vertices[v2].Q[2] += area * (n.x * n.z);
                vertices[v2].Q[3] += area * (n.x * d);
                vertices[v2].Q[4] += area * (n.y * n.y);
                vertices[v2].Q[5] += area * (n.y * n.z);
                vertices[v2].Q[6] += area * (n.y * d);
                vertices[v2].Q[7] += area * (n.z * n.z);
                vertices[v2].Q[8] += area * (n.z * d);
                vertices[v2].Q[9] += area * (d * d);
            }
        }
    }

    private func manifoldConstraint(v1: Int, v2: Int) -> Bool {
        var vertices = Set<Int>()

        var idEdgeV1V2 = 0
        for idEdge1 in self.vertices[v1].edges {
            let a = (edges[idEdge1].v1 == v1) ? edges[idEdge1].v2 : edges[idEdge1].v1
            vertices.insert(a)
            if (a != v2)
            {
                for idEdge2 in self.vertices[v2].edges {
                    let b = (edges[idEdge2].v1 == v2) ? edges[idEdge2].v2 : edges[idEdge2].v1;
                    vertices.insert(b)
                    if a == b && getTriangle(v1: v1, v2: v2, v3: a) == -1 {
                        return false;
                    }
                }
            } else {
                idEdgeV1V2 = idEdge1;
            }
        }
        if vertices.count <= 4 || ( self.vertices[v1].onBoundary && self.vertices[v2].onBoundary && !edges[idEdgeV1V2].onBoundary) {
            return false;
        }
        return true;
    }



    private func computeEdgeCost(v1: Int, v2: Int, newPos: inout Vector3) -> Double {
        var Q = [Double](repeating: 0, count: 10)
        var M = [Double](repeating: 0, count: 12)
        var pos = Vector3(x: 0, y: 0, z: 0)
        for i in 0..<10{
            Q[i] = vertices[v1].Q[i] + vertices[v2].Q[i];
        }
        M[0] = Q[0]; // (0, 0)
        M[1] = Q[1]; // (0, 1)
        M[2] = Q[2]; // (0, 2)
        M[3] = Q[3]; // (0, 3)
        M[4] = Q[1]; // (1, 0)
        M[5] = Q[4]; // (1, 1)
        M[6] = Q[5]; // (1, 2)
        M[7] = Q[6]; // (1, 3)
        M[8] = Q[2]; // (2, 0)
        M[9] = Q[5]; // (2, 1)
        M[10] = Q[7]; // (2, 2);
        M[11] = Q[8]; // (2, 3);
        let det =   M[0] * M[5] * M[10] + M[1] * M[6] * M[8] + M[2] * M[4] * M[9] - M[0] * M[6] * M[9]  - M[1] * M[4] * M[10] - M[2] * M[5] * M[8];
        if det != 0.0 {
            let d = 1.0 / det;
            pos.x = d * (M[1]*M[7]*M[10] + M[2]*M[5]*M[11] + M[3]*M[6]*M[9] - M[1]*M[6]*M[11] - M[2]*M[7]*M[9]  - M[3]*M[5]*M[10]);
            pos.y = d * (M[0]*M[6]*M[11] + M[2]*M[7]*M[8]  + M[3]*M[4]*M[10] - M[0]*M[7]*M[10] - M[2]*M[4]*M[11] - M[3]*M[6]*M[8]);
            pos.z = d * (M[0]*M[7]*M[9]  + M[1]*M[4]*M[11] + M[3]*M[5]*M[8]  - M[0]*M[5]*M[11] - M[1]*M[7]*M[8]  - M[3]*M[4]*M[9]);
            newPos.x = pos.x
            newPos.y = pos.y
            newPos.z = pos.z
        } else {
            let w = 0.5
            newPos = w * points[v1] + w * points[v2];
            pos.x = newPos.x
            pos.y = newPos.y
            pos.z = newPos.z
        }

        let qem = pos.x  * (Q[0] * pos.x + Q[1] * pos.y + Q[2] * pos.z + Q[3]) +
            pos.y  * (Q[1] * pos.x + Q[4] * pos.y + Q[5] * pos.z + Q[6]) +
            pos.z  * (Q[2] * pos.x + Q[5] * pos.y + Q[7] * pos.z + Q[8]) +
            (Q[3] * pos.x + Q[6] * pos.y + Q[8] * pos.z + Q[9])

        let oldPosV1 =  points[v1];
        let oldPosV2 =  points[v2];
        var tris = vertices[v1].triangles
        for idTriangle in vertices[v2].triangles {
            tris.append(idTriangle);
        }

        for idTriangle in tris {
            let a0 = triangles[idTriangle].v0
            let a1 = triangles[idTriangle].v1
            let a2 = triangles[idTriangle].v2

            let n1 = ((points[a1] - points[a0]) ^ (points[a2] - points[a0])).normalized

            points[v1] = newPos;
            points[v2] = newPos;

            let n2 = (points[a1] - points[a0]) ^ (points[a2] - points[a0]).normalized

            points[v1] = oldPosV1;
            points[v2] = oldPosV2;

            if n1 * n2 < 0.0 {
                return Double.greatestFiniteMagnitude
            }
        }

        if ecolManifoldConstraint && !manifoldConstraint(v1: v1, v2: v2) {
            return Double.greatestFiniteMagnitude
        }
        return qem;

    }

    private func edgeCollapse(qem: inout Double) -> Bool {
        var currentEdge = EdgePriorityQueue()

        var done = false;
        repeat {
            done = false;
            if pqueue.isEmpty {
                done = true;
                break;
            } else {
                currentEdge = pqueue.pop()!

            }
        } while ( (!edges[currentEdge.name].tag) || (edges[currentEdge.name].qem != currentEdge.qem))

        if done {
            return false;
        }

        let v1 = edges[currentEdge.name].v1;
        let v2 = edges[currentEdge.name].v2;

        qem = currentEdge.qem;
        edgeCollapse(v1: v1, v2: v2);
        points[v1] = edges[currentEdge.name].position
        for k in 0..<10 {
            vertices[v1].Q[k] += vertices[v2].Q[k];
        }

        // Update priority queue
        var incidentVertices = [Int]()
        for idEdge in vertices[v1].edges {
            let a = edges[idEdge].v1;
            let b = edges[idEdge].v2;
            incidentVertices.append((a != v1) ? a : b)
            var pqEdge = EdgePriorityQueue(name: idEdge, qem: 0)
            let qem =  computeEdgeCost(v1: a, v2: b, newPos: &edges[idEdge].position);
            edges[idEdge].qem = qem
            pqEdge.qem = qem
            pqueue.push(pqEdge);
        }

        for idVertex in incidentVertices {
            for idEdge in vertices[idVertex].edges {
                let a = edges[idEdge].v1;
                let b = edges[idEdge].v2;
                if a != v1 && b != v1 {
                    var pqEdge = EdgePriorityQueue(name: idEdge, qem: 0)
                    let qem =  computeEdgeCost(v1: a, v2: b, newPos: &edges[idEdge].position);
                    edges[idEdge].qem = qem
                    pqEdge.qem = qem
                    pqueue.push(pqEdge);
                }
            }
        }
        return true;
    }


    func decimate( targetVertices: Int, targetTriangles: Int, targetError: Double) {
        var qem = 0.0;

        initializeQEM();
        initializePriorityQueue();

        let invDiag = 1 / diagBB;
        while !pqueue.isEmpty && edgeCount > 0 && vertCount > targetVertices && triangleCount > targetTriangles && qem < targetError {

            if !edgeCollapse(qem: &qem) {
                break;
            }
            if qem < 0.0 {
                qem = 0.0;
            } else {
                qem = sqrt(qem) * invDiag;
            }
        }

    }

    func getMeshData() -> ([Vector3], [Triangle]) {
        var map = [Int](repeating: 0, count: pointCount)
        var points = [Vector3]()
        var triangles = [Triangle]()
        var counter = 0;
        for v in 0..<pointCount {
            if vertices[v].tag {
                points.append(self.points[v])
                map[v] = counter
                counter += 1
            }
        }
        counter = 0;
        for t in 0..<initialTriangles {
            if trianglesTags[t] {
                let triangle = Triangle(v0:  map[self.triangles[t].v0], v1: map[self.triangles[t].v1], v2: map[self.triangles[t].v2])
                triangles.append(triangle)
            }
        }
        return (points, triangles)
    }


    struct Vertex {
        var edges = [Int]()
        var triangles = [Int]()
        var Q: [Double] = [0,0,0,0,0,0,0,0,0,0,0]
        // 0 1 2 3
        //   4 5 6
        //     7 8
        //       9
        var tag = true
        var onBoundary = false

        mutating func remove(triangle: Int) {
            if let index = triangles.index(of: triangle) {
                triangles.remove(at: index)
            }
        }

        mutating func remove(edge: Int) {
            if let index = edges.index(of: edge) {
                edges.remove(at: index)
            }
        }
        mutating func resetQ() {
            memset(&Q, 0, Q.count * MemoryLayout<Double>.size)
        }
    }

    struct Edge {
        var v1 = 0
        var v2 = 0
        var qem = 0.0
        var position = Vector3(x: 0, y: 0,z: 0)
        var tag = false
        var onBoundary = false
    }
    struct EdgePriorityQueue: Comparable {
        var name = 0
        var qem = 0.0
    }
}
func ==(lhs:Decimator.EdgePriorityQueue, rhs: Decimator.EdgePriorityQueue) -> Bool { return (lhs.name == rhs.name ) }

func <(lhs:Decimator.EdgePriorityQueue, rhs: Decimator.EdgePriorityQueue) -> Bool { return (lhs.qem > rhs.qem) }
func >(lhs:Decimator.EdgePriorityQueue, rhs: Decimator.EdgePriorityQueue) -> Bool { return (lhs.qem < rhs.qem) }

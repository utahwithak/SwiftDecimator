//
//  Decimator.swift
//  Swift Decimator
//
//  Created by Carl Wieland on 1/24/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation

public struct Triangle {
    var v0 = 0
    var v1 = 0
    var v2 = 0
}

public class Decimator {

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

    public convenience init?(OBJAt file: URL) {
        var points = [Vector3]()
        var triangles = [Triangle]()


        guard let text = try? String(contentsOf: file) else {
            return nil
        }
        let lines = text.components(separatedBy: "\n")

        for buffer in lines {
            let components = buffer.components(separatedBy: " ")
            if components.first == "v" {
                guard let x = Double(components[1]),
                    let y =  Double(components[2]),
                    let z = Double(components[3]) else {
                        print("invalid V")
                        return nil
                }
                points.append(Vector3(x,y,z))
            } else if components.first == "f" {
                guard let v0 = Int(components[1]),
                    let v1 =  Int(components[2]),
                    let v2 = Int(components[3]) else {
                        print("invalid F")
                        return nil
                }

                triangles.append(Triangle(v0: v0 - 1, v1: v1 - 1, v2: v2 - 1))
            }

        }

        guard !points.isEmpty && !triangles.isEmpty else {
            return nil
        }

        self.init(points: points, triangles: triangles)
    }

    public init( points: [Vector3], triangles: [Triangle]) {
        self.points = points
        self.triangles = triangles
        initialTriangles = triangles.count
        triangleCount = triangles.count
        vertCount = points.count
        edges = [Edge]()
        edges.reserveCapacity(triangleCount * 3)
        vertices = [Vertex](repeating:Vertex(), count: vertCount)

        var edge = Edge()
        edge.tag = true;
        edge.onBoundary = true;

        for t in 0..<triangleCount {
            let tri = [triangles[t].v0, triangles[t].v1, triangles[t].v2]
            trianglesTags.append(true)
            for k in 0..<3 {
                edge.v1 = tri[k];
                edge.v2 = tri[(k+1) % 3];
                if !vertices[edge.v1].triangles.contains(t) {
                    vertices[edge.v1].triangles.append(t)
                }
                if let idEdge = getEdge(v1: edge.v1, v2: edge.v2) {
                    edges[idEdge].onBoundary = false;
                } else {
                    edges.append(edge);
                    if !vertices[edge.v1].edges.contains(edgeCount) {
                        vertices[edge.v1].edges.append(edgeCount);
                    }
                    if !vertices[edge.v2].edges.contains(edgeCount) {
                        vertices[edge.v2].edges.append(edgeCount);
                    }
                    edgeCount += 1
                }
            }
        }

        edgeCount = edges.count

        for v in 0..<vertCount {
            vertices[v].onBoundary = vertices[v].edges.contains(where: { edges[$0].onBoundary })
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
            return (edges[idEdge].v1 == v1 && edges[idEdge].v2 == v2) || (edges[idEdge].v1 == v2 && edges[idEdge].v2 == v1)
        })
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
            } else if getTriangle(v1: v1, v2: u, v3: w) == -1 {

                vertices[v1].insert(triangle: idTriangle)

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

        var it = 0
        var idEdge = 0
        while it < vertices[v2].edges.count {

            idEdge = vertices[v2].edges[it]
            let w = (edges[idEdge].v1 == v2) ? edges[idEdge].v2 : edges[idEdge].v1;
            if w == v1 {
                edges[idEdge].tag = false;
                vertices[w].remove(edge: idEdge)
                edgeCount -= 1
            } else if getEdge(v1: v1, v2: w) == nil {
                if edges[idEdge].v1 == v2 {
                    edges[idEdge].v1 = v1;
                } else{
                    edges[idEdge].v2 = v1;
                }

                vertices[v1].insert(edge: idEdge)
            } else {
                edges[idEdge].tag = false;
                vertices[w].remove(edge: idEdge)
                edgeCount -= 1
            }

            it += 1
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
            vertices[idVertex].onBoundary = vertices[idVertex].edges.contains(where: { edges[$0].onBoundary})
        }
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

        for v in 0..<pointCount {
            for idTriangle in vertices[v].triangles {
                let i = triangles[idTriangle].v0
                let j = triangles[idTriangle].v1
                let k = triangles[idTriangle].v2
                var n = (points[j] - points[i]) ^ (points[k] - points[i]);
                let area = n.magnitude

                n.normalize()
                let d = -(points[v] * n)
                vertices[v].Q0 += area * (n.x * n.x);
                vertices[v].Q1 += area * (n.x * n.y);
                vertices[v].Q2 += area * (n.x * n.z);
                vertices[v].Q3 += area * (n.x * d);
                vertices[v].Q4 += area * (n.y * n.y);
                vertices[v].Q5 += area * (n.y * n.z);
                vertices[v].Q6 += area * (n.y * d);
                vertices[v].Q7 += area * (n.z * n.z);
                vertices[v].Q8 += area * (n.z * d);
                vertices[v].Q9 += area * (d     * d);
            }
        }

        let w = 1000.0
        for edge in edges {
            let v1 = edge.v1
            let v2 = edge.v2
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
                let area = w * (u1^u2).magnitude
                u1.normalize()
                let n = (u2 - (u2 * u1) * u1).normalized


                var d = -(points[v1] * n);
                vertices[v1].Q0 += area * (n.x * n.x);
                vertices[v1].Q1 += area * (n.x * n.y);
                vertices[v1].Q2 += area * (n.x * n.z);
                vertices[v1].Q3 += area * (n.x * d);
                vertices[v1].Q4 += area * (n.y * n.y);
                vertices[v1].Q5 += area * (n.y * n.z);
                vertices[v1].Q6 += area * (n.y * d);
                vertices[v1].Q7 += area * (n.z * n.z);
                vertices[v1].Q8 += area * (n.z * d);
                vertices[v1].Q9 += area * (d * d);

                d = -(points[v2] * n);
                vertices[v2].Q0 += area * (n.x * n.x);
                vertices[v2].Q1 += area * (n.x * n.y);
                vertices[v2].Q2 += area * (n.x * n.z);
                vertices[v2].Q3 += area * (n.x * d);
                vertices[v2].Q4 += area * (n.y * n.y);
                vertices[v2].Q5 += area * (n.y * n.z);
                vertices[v2].Q6 += area * (n.y * d);
                vertices[v2].Q7 += area * (n.z * n.z);
                vertices[v2].Q8 += area * (n.z * d);
                vertices[v2].Q9 += area * (d * d);
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

        let Q0 = vertices[v1].Q0 + vertices[v2].Q0
        let Q1 = vertices[v1].Q1 + vertices[v2].Q1
        let Q2 = vertices[v1].Q2 + vertices[v2].Q2
        let Q3 = vertices[v1].Q3 + vertices[v2].Q3
        let Q4 = vertices[v1].Q4 + vertices[v2].Q4
        let Q5 = vertices[v1].Q5 + vertices[v2].Q5
        let Q6 = vertices[v1].Q6 + vertices[v2].Q6
        let Q7 = vertices[v1].Q7 + vertices[v2].Q7
        let Q8 = vertices[v1].Q8 + vertices[v2].Q8
        let Q9 = vertices[v1].Q9 + vertices[v2].Q9


        let q22 = Q2 * Q2
        let q11 = Q1 * Q1
        let q55 = Q5 * Q5
        let q15 = Q1 * Q5
        let q34 = Q3 * Q4
        let q06 = Q0 * Q6
        let det = Q0 * Q4 * Q7 +
            q15 * Q2 +
            Q2 * Q1 * Q5 -
            Q0 * q55 -
            q11 * Q7 -
            q22 * Q4 ;

        var pos = Vector3(x: 0, y: 0, z: 0)

        if det != 0.0 {
            let d = 1.0 / det;
            pos.x = d * (Q1 * Q6 * Q7 + Q2 * Q4 * Q8 + Q3 * q55 - q15 * Q8 - Q2 * Q6 * Q5  - q34 * Q7);
            pos.y = d * (Q0 * Q5 * Q8 + q22 * Q6  + Q3 * Q1 * Q7 - q06 * Q7 - Q2 * Q1 * Q8 - Q3 * Q5 * Q2);
            pos.z = d * (q06 * Q5  + q11 * Q8 + q34 * Q2  - Q0 * Q4 * Q8 - Q1 * Q6 * Q2  - Q3 * q15);
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

        let qem = pos.x  * (Q0 * pos.x + Q1 * pos.y + Q2 * pos.z + Q3) +
            pos.y  * (Q1 * pos.x + Q4 * pos.y + Q5 * pos.z + Q6) +
            pos.z  * (Q2 * pos.x + Q5 * pos.y + Q7 * pos.z + Q8) +
            (Q3 * pos.x + Q6 * pos.y + Q8 * pos.z + Q9)

        let oldPosV1 =  points[v1];
        let oldPosV2 =  points[v2];




        for idTriangle in vertices[v1].triangles {
            let a0 = triangles[idTriangle].v0
            let a1 = triangles[idTriangle].v1
            let a2 = triangles[idTriangle].v2

            let n1 = ((points[a1] - points[a0]) ^ (points[a2] - points[a0]) )//.normalized

            points[v1] = newPos;
            points[v2] = newPos;

            let n2 =  ((points[a1] - points[a0]) ^ (points[a2] - points[a0]))//.normalized

            points[v1] = oldPosV1;
            points[v2] = oldPosV2;

            if n1 * n2 < 0.0 {
                return Double.greatestFiniteMagnitude
            }
        }
        for idTriangle in vertices[v2].triangles where !vertices[v1].triangles.contains(idTriangle) {
            let a0 = triangles[idTriangle].v0
            let a1 = triangles[idTriangle].v1
            let a2 = triangles[idTriangle].v2

            let n1 = ((points[a1] - points[a0]) ^ (points[a2] - points[a0]) )//.normalized

            points[v1] = newPos;
            points[v2] = newPos;

            let n2 =  ((points[a1] - points[a0]) ^ (points[a2] - points[a0]))//.normalized

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

        vertices[v1].Q0 += vertices[v2].Q0
        vertices[v1].Q1 += vertices[v2].Q1
        vertices[v1].Q2 += vertices[v2].Q2
        vertices[v1].Q3 += vertices[v2].Q3
        vertices[v1].Q4 += vertices[v2].Q4
        vertices[v1].Q5 += vertices[v2].Q5
        vertices[v1].Q6 += vertices[v2].Q6
        vertices[v1].Q7 += vertices[v2].Q7
        vertices[v1].Q8 += vertices[v2].Q8
        vertices[v1].Q9 += vertices[v2].Q9


        // Update priority queue
        var incidentVertices = [Int]()
        for idEdge in vertices[v1].edges {
            let a = edges[idEdge].v1;
            let b = edges[idEdge].v2;
            incidentVertices.append((a != v1) ? a : b)
            var pqEdge = EdgePriorityQueue(name: idEdge, qem: 0)
            let qem = computeEdgeCost(v1: a, v2: b, newPos: &edges[idEdge].position);
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
                    let qem = computeEdgeCost(v1: a, v2: b, newPos: &edges[idEdge].position);
                    edges[idEdge].qem = qem
                    pqEdge.qem = qem
                    pqueue.push(pqEdge);
                }
            }
        }
        return true;
    }


    public func decimate( targetVertices: Int, targetTriangles: Int, targetError: Double) {
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

    public func getMeshData() -> ([Vector3], [Triangle]) {
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


        var Q0 = 0.0
        var Q1 = 0.0
        var Q2 = 0.0
        var Q3 = 0.0
        var Q4 = 0.0
        var Q5 = 0.0
        var Q6 = 0.0
        var Q7 = 0.0
        var Q8 = 0.0
        var Q9 = 0.0

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
        mutating func insert(triangle: Int) {
            guard !triangles.contains(triangle) else {
                return
            }
            triangles.append(triangle)
        }
        mutating func insert(edge: Int) {
            guard !edges.contains(edge) else {
                return
            }
            edges.append(edge)
        }
    }

    public func generateOBJ() -> String {
        var file = ""
        let (points,triangles) = getMeshData()

        for p in points {
            file += "v \(p.x) \(p.y) \(p.z)\n"
        }
        
        for tri in triangles {
            file += "f \(tri.v0 + 1) \(tri.v1 + 1) \(tri.v2 + 1)\n"
        }

        return file
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

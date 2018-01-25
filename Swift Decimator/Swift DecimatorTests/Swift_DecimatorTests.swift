//
//  Swift_DecimatorTests.swift
//  Swift DecimatorTests
//
//  Created by Carl Wieland on 1/24/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import XCTest
import simd
@testable import Swift_Decimator

class Swift_DecimatorTests: XCTestCase {
    
    override func setUp() {
        super.setUp()
        // Put setup code here. This method is called before the invocation of each test method in the class.
    }
    
    override func tearDown() {
        // Put teardown code here. This method is called after the invocation of each test method in the class.
        super.tearDown()
    }
    
    func testExample() {
        let (positions, faces) = IcoSphereCreator.icospherePoints(level: 3, radius: 10, shift: float3(0,0,0))

        let rawPoints = positions.map{ Vector3($0)}
        let decimator = Decimator(points: rawPoints, triangles: faces)
        decimator.decimate(targetVertices: 300, targetTriangles: 600, targetError: 1.0)
        let (points, tris) = decimator.getMeshData()
        print("points:\(points.count) tris:\(tris.count)")
    }
    
    func testPerformanceExample() {
        // This is an example of a performance test case.
        self.measure {
            // Put the code you want to measure the time of here.
        }
    }
    
}

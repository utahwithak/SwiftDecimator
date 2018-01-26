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
    
//    func testExample() {
//        let (positions, faces) = IcoSphereCreator.icospherePoints(level: 3, radius: 10, shift: float3(0,0,0))
//
//        let rawPoints = positions.map{ Vector3($0)}
//        let decimator = Decimator(points: rawPoints, triangles: faces)
//        decimator.decimate(targetVertices: 300, targetTriangles: 600, targetError: 1.0)
//        let (points, tris) = decimator.getMeshData()
//        print("points:\(points.count) tris:\(tris.count)")
//    }

    func testOBJ() {
        let bundle = Bundle(for: type(of:self))
        guard let url = bundle.url(forResource: "pumpkin_tall_10k", withExtension: "obj") else {
            XCTAssert(false)
            return
        }

        guard let decimator = Decimator(OBJAt: url) else {
            XCTAssert(false, "invalid url")
            return
        }

        decimator.decimate(targetVertices: 100, targetTriangles: 200, targetError: 1)
        let objFile = decimator.generateOBJ()
        guard !objFile.isEmpty else {
            XCTAssert(false)
            return
        }
//
//        // The template string:
//        let template = URL(fileURLWithPath: NSTemporaryDirectory()).appendingPathComponent("out.XXXXXX") as NSURL
//
//        // Fill buffer with a C string representing the local file system path.
//        var buffer = [Int8](repeating: 0, count: Int(PATH_MAX))
//        template.getFileSystemRepresentation(&buffer, maxLength: buffer.count)
//
//        // Create unique file name (and open file):
//        let fd = mkstemp(&buffer)
//        if fd != -1 {
//
//            // Create URL from file system string:
//            let url = NSURL(fileURLWithFileSystemRepresentation: buffer, isDirectory: false, relativeTo: nil).appendingPathExtension("obj")!
//            print(url.path)

            let path = "/Users/carl8382/Desktop/swift.obj"

        try! objFile.write(to: URL(fileURLWithPath: path), atomically: true, encoding: .utf8)
//
//        } else {
//            print("Error: " + String(cString: strerror(errno)))
//        }
    }
    
    func testPerformanceExample() {
        // This is an example of a performance test case.
        self.measure {
            // Put the code you want to measure the time of here.
        }
    }
    
}

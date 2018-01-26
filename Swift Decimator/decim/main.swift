//
//  main.swift
//  decim
//
//  Created by Carl Wieland on 1/26/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation

print("Hello, World!")

guard CommandLine.arguments.count == 6 else {
    print( "Usage: ./MeshSimplification fileName targetNTrianglesDecimatedMesh targetNVerticesDecimatedMesh maxDecimationError fileNameOut.obj")
        exit(EXIT_FAILURE)
}

let argv = CommandLine.arguments
let fileName = argv[1]
guard let targetNTrianglesDecimatedMesh = Int(argv[2]),
      let targetNVerticesDecimatedMesh  = Int(argv[3]),
      let maxDecimationError = Double(argv[4]) else {
        print("Invalid arguments")
        exit(EXIT_FAILURE)

}
let fileNameOut = argv[5]

for _ in 0..<10 {

    guard let decimator = Decimator(OBJAt: URL(fileURLWithPath: fileName)) else {
        print("Failed to create OBJ")
        exit(EXIT_FAILURE)

    }
    decimator.decimate(targetVertices: targetNVerticesDecimatedMesh, targetTriangles: targetNTrianglesDecimatedMesh, targetError: maxDecimationError)

    let objText = decimator.generateOBJ()

    do {
        try objText.write(to: URL(fileURLWithPath: fileNameOut), atomically: true, encoding: .utf8)
    } catch {
        print("error")
        exit(EXIT_FAILURE)
    }
}
exit(EXIT_SUCCESS)

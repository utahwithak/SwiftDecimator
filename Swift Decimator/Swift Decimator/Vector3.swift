//
//  Vector3.swift
//  Swift Decimator
//
//  Created by Carl Wieland on 1/24/18.
//  Copyright © 2018 Datum Apps. All rights reserved.
//

import Foundation
import simd

public struct Vector3 {
    public static let zero = Vector3(x: 0, y: 0, z: 0)

    var x: Double
    var y: Double
    var z: Double

    init(){
        x = 0
        y = 0
        z = 0
    }
    init(_ x: Double, _ y: Double, _ z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
    init(_ val: float3) {
        self.x = Double(val.x)
        self.y = Double(val.y)
        self.z = Double(val.z)
    }
    init( x: Double,  y: Double,  z: Double) {
        self.x = x
        self.y = y
        self.z = z
    }
    var magnitude: Double {
        return sqrt(magnitudeSquared)
    }

    var magnitudeSquared: Double {
        return (x * x) + (y * y) + (z * z)
    }

    mutating func normalize() {
        if magnitude != 0 {
            self /= magnitude
        }
    }
    var normalized: Vector3 {
        let magnitude = self.magnitude
        if magnitude != 0 {
            return self / magnitude
        }
        return self
    }
    func dot(_ other: Vector3) -> Double {
        return self * other
    }
    var asFloat3: float3 {
        return float3(Float(x), Float(y), Float(z))
    }
}

func /=(lhs: inout Vector3, rhs: Double) {
    lhs.x /= rhs
    lhs.y /= rhs
    lhs.z /= rhs
}

func *=(lhs: inout Vector3, rhs: Double) {
    lhs.x *= rhs
    lhs.y *= rhs
    lhs.z *= rhs
}

func -=(lhs: inout Vector3, rhs: Vector3) {
    lhs.x -= rhs.x
    lhs.y -= rhs.y
    lhs.z -= rhs.z
}
func +=(lhs: inout Vector3, rhs: Vector3) {
    lhs.x += rhs.x
    lhs.y += rhs.y
    lhs.z += rhs.z
}

func +(lhs: Vector3, rhs: Vector3) -> Vector3 {
    return Vector3( x: lhs.x + rhs.x, y: lhs.y + rhs.y,z: lhs.z + rhs.z)
}

prefix func -(lhs:inout Vector3){
    lhs.x *= -1
    lhs.y *= -1
    lhs.z *= -1
}

prefix func -(lhs: Vector3) -> Vector3{
    return Vector3(x: lhs.x * -1, y: lhs.y * -1, z: lhs.z * -1)
}
func -(lhs: Vector3, rhs: Vector3) -> Vector3 {
    return Vector3( x: lhs.x - rhs.x, y: lhs.y - rhs.y,z: lhs.z - rhs.z)
}

func *(lhs: Vector3, rhs: Double) -> Vector3 {
    return Vector3( x: lhs.x * rhs, y: lhs.y * rhs,z: lhs.z * rhs)
}

func *(rhs: Double, lhs: Vector3) -> Vector3 {
    return Vector3( x: lhs.x * rhs, y: lhs.y * rhs,z: lhs.z * rhs)
}

func /(lhs: Vector3, rhs: Double) -> Vector3 {
    return Vector3( x: lhs.x / rhs, y: lhs.y / rhs,z: lhs.z / rhs)
}

func *(lhs: Vector3, rhs: Vector3) -> Double {
    return lhs.x * rhs.x + lhs.y * rhs.y + lhs.z * rhs.z
}

func ^(lhs: Vector3, rhs: Vector3) -> Vector3 {

    return Vector3(x: lhs.y * rhs.z - lhs.z * rhs.y,
                   y: lhs.z * rhs.x - lhs.x * rhs.z,
                   z: lhs.x * rhs.y - lhs.y * rhs.x)
}


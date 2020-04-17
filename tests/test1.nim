import math, random, unittest
import kdtree

suite "kdtree test suite":
  
  let numPoints = 10_000

  proc getTree(): KdTree[2, int] = 
    var
      pointData = newSeqOfCap[(array[2, float], int)](numPoints)
      x: float
      y: float
      r = initRand(34)

    for a in 0..<numPoints:
      x = r.rand(100.0)
      y = r.rand(100.0)
      pointData.add(([x, y], a))

    result = newKdTree[2, int](pointData)

  proc getSimpleTree(): KdTree[2, int] =
    let pointsData = [([2.0, 3.0], 1), 
                          ([5.0, 4.0], 2), 
                          ([9.0, 6.0], 3), 
                          ([4.0, 7.0], 4), 
                          ([8.0, 1.0], 5), 
                          ([7.0, 2.0], 6)]

    result = newKdTree[2, int](pointsData)

  test "Build tree":
    var
      points = newSeqOfCap[array[2, float]](numPoints)
      values = newSeqOfCap[int](numPoints)
      x: float
      y: float
      r = initRand(34)

    for a in 0..<numPoints:
      x = r.rand(100.0)
      y = r.rand(100.0)
      points.add([x, y])
      values.add(a)

    var tree = newKdTree[2, int](points, values)

    check(tree.len() == numPoints)
    check(tree.height() == 14)
    
    expect AssertionError:
      discard points.pop()
      discard newKdTree[2, int](points, values)

    expect AssertionError:
      points.setLen(0)
      discard newKdTree[2, int](points, values)

  test "Balance tree":
    var 
      tree = getTree()
      numPoints = tree.len()
      r = initRand(34)
      x: float
      y: float
      
    for a in 0..<1000:
      x = r.rand(20.0)
      y = r.rand(20.0)
      tree.add([x, y], numPoints+a)

    check(tree.isBalanced() == 8)

    tree.rebalance()
    check(tree.isBalanced == 0)

  test "Custom distance function":
    var
      points = newSeqOfCap[array[2, float]](numPoints)
      values = newSeqOfCap[int](numPoints)
      x: float
      y: float
      r = initRand(34)

    proc myDistFunc(p1, p2: KdPoint[2]): float {.closure.} =
      result = 0.0
      for i in 0..<len(p1):
        result += (p1[i] - p2[i]) * (p1[i] - p2[i])
      result = sqrt(result)

    for a in 0..<numPoints:
      x = r.rand(100.0)
      y = r.rand(100.0)
      points.add([x, y])
      values.add(a)

    var tree = newKdTree[2, int](points, values, distFunc=myDistFunc)

    let (_, value, dist) = tree.nearestNeighbour([50.0, 50.0])
    check:
      value == 922 
      abs(dist - 0.4433013003347581) <= 1E-7

  test "Nearest-neighbour search":
    var tree = getTree()
    let (_, value, dist) = tree.nearestNeighbour([50.0, 50.0])
    check:
      value == 922 
      abs(dist - 0.1965160428784874) <= 1E-7

    var tree2 = getSimpleTree()
    let (_, value2, dist2) = tree2.nearestNeighbour([9.0, 2.0])
    check:
      value2 == 5 
      abs(dist2 - 2.0) <= 1E-7

  test "Nearest k neighbours search":
    var tree = getTree()
    let ret = tree.nearestNeighbours([50.0, 50.0], numNeighbours=5)

    check:
      ret[0][1] == 922
      ret[1][1] == 3110
      ret[2][1] == 2952
      ret[3][1] == 8600
      ret[4][1] == 2899

      abs(ret[0][2] - 0.1965160428784874) <= 1E-7
      abs(ret[1][2] - 0.7045863318354959) <= 1E-7
      abs(ret[2][2] - 1.131211308200418) <= 1E-7
      abs(ret[3][2] - 2.01025020458314) <= 1E-7
      abs(ret[4][2] - 2.27941022920038) <= 1E-7

    # for (pt, value, dist) in ret:
    #   echo pt, ", ", value, ", ", dist

    var tree2 = getSimpleTree()
    let ret2 = tree2.nearestNeighbours([9.0, 2.0], 3)
    check:
      ret2[0][1] == 5
      ret2[1][1] == 6
      ret2[2][1] == 3

      abs(ret2[0][2] - 2.0) <= 1E-7
      abs(ret2[1][2] - 4.0) <= 1E-7
      abs(ret2[2][2] - 16.0) <= 1E-7

    # for (pt, value, dist) in ret2:
    #   echo pt, ", ", value, ", ", dist

  test "Within-radius search":
    var tree = getTree()
    let ret = tree.withinRadius([50.0, 50.0], radius=1.5, sortResults=true)
    
    check:
      ret[0][1] == 922
      ret[1][1] == 3110
      ret[2][1] == 2952
      ret[3][1] == 8600

      abs(ret[0][2] - 0.1965160428784874) <= 1E-7
      abs(ret[1][2] - 0.7045863318354959) <= 1E-7
      abs(ret[2][2] - 1.131211308200418) <= 1E-7
      abs(ret[3][2] - 2.01025020458314) <= 1E-7

    # for (pt, value, dist) in ret:
    #   echo pt, ", ", value, ", ", dist

    var tree2 = getSimpleTree()
    let ret2 = tree2.withinRadius([9.0, 2.0], radius=3.0, sortResults=true)
    check:
      ret2[0][1] == 5
      ret2[1][1] == 6

      abs(ret2[0][2] - 2.0) <= 1E-7
      abs(ret2[1][2] - 4.0) <= 1E-7

    # for (pt, value, dist) in ret2:
    #   echo pt, ", ", value, ", ", dist

  test "Within-range test":
    var
      tree = getTree()
      min = [50.0-1.5, 50.0-1.5]
      max = [50.0+1.5, 50.0+1.5]
      hyperRect = newHyperRectangle(min, max)

    var ret = tree.withinRange(hyperRect)
    # for (pt, value) in ret:
    #   echo "point=", pt, ", value=", value

    check:
      ret[0][1] == 260
      ret[1][1] == 9842
      ret[2][1] == 2952
      ret[3][1] == 8600
      ret[4][1] == 3110
      ret[5][1] == 2899
      ret[6][1] == 922

    var tree2 = getSimpleTree()
    min = [4.0, 1.0]
    max = [9.0, 5.0]
    hyperRect = newHyperRectangle(min, max)
    let ret2 = tree2.withinRange(hyperRect)
    # for (pt, value) in ret2:
    #   echo "point=", pt, ", value=", value
    check:
      ret2[0][1] == 6
      ret2[1][1] == 5
      ret2[2][1] == 2

  test "3D tree":
    var
      pointData2 = [ 
        ([1.0, 1.0, 1.1], 1), 
        ([1.0, -1.0, 1.2], 2), 
        ([-1.0, -1.0, 1.3], 3), 
        ([-1.0, 1.0, 1.4], 4), 
        ([1.0, 1.0, -1.4], 5), 
        ([1.0, -1.0, -1.3], 6), 
        ([-1.0, -1.0, -1.2], 7), 
        ([-1.0, 1.0, -1.1], 8), 
        ([0.5, -0.5, 0.5], 9)]
      

    var tree3D = newKdTree[3, int](pointData2)

    let (_, value, dist) = tree3D.nearestNeighbour([0.0, 0.0, 0.0])
    # echo "result"
    # echo pt
    # echo value
    # echo dist
    check:
      value == 9 
      abs(dist - 0.75) <= 1E-7

    var ret = tree3D.nearestNeighbours([0.0, 0.0, 0.0], 9)
    # for (pt, value, dist) in ret:
    #   echo pt
    #   echo value
    #   echo dist
    
    check:
      ret[0][1] == 9 
      ret[1][1] == 8
      ret[2][1] == 1 
      ret[3][1] == 7 
      ret[4][1] == 2 
      ret[5][1] == 3 
      ret[6][1] == 6 
      ret[7][1] == 4 
      ret[8][1] == 5
      abs(ret[0][2] - 0.75) <= 1E-7
      abs(ret[1][2] - 3.21) <= 1E-7
      abs(ret[2][2] - 3.21) <= 1E-7
      abs(ret[3][2] - 3.44) <= 1E-7
      abs(ret[4][2] - 3.44) <= 1E-7
      abs(ret[5][2] - 3.69) <= 1E-7
      abs(ret[6][2] - 3.69) <= 1E-7
      abs(ret[7][2] - 3.96) <= 1E-7
      abs(ret[8][2] - 3.96) <= 1E-7

  
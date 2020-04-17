## kdtree is a pure Nim k-d tree implementation. k-d trees are data structures for performing
## efficient spatial query operations on point data sets.
## 
## .. code-block:: nim
##   import random, strformat
##   import kdtree
## 
##   let numPoints = 100_000
##   var
##     points = newSeqOfCap[array[2, float]](numPoints)
##     values = newSeqOfCap[int](numPoints)
##     x: float
##     y: float
##     r = initRand(34)
## 
##   for a in 0..<numPoints:
##     x = r.rand(100.0)
##     y = r.rand(100.0)
##     points.add([x, y])
##     values.add(a)
## 
##   echo fmt"Building tree of {numPoints} random points..."
## 
##   # Create a 2-D tree of int's
##   var tree = newKdTree[2, int](points, values)
## 
##   # Perform nearestNeighour searches
##   let numSearches = 10_000
##   for a in 0..<numSearches:
##     x = r.rand(100.0)
##     y = r.rand(100.0)
##     let (pt, values, dist) = tree.nearestNeighbour([x, y])
##     echo fmt"point={pt}, value={value}, dist={dist}"
## 
##   # Perform nearestNeighours searches
##   let n = 10
##   for a in 0..<numSearches:
##     x = r.rand(100.0)
##     y = r.rand(100.0)
##     let ret = tree.nearestNeighbours([x, y], n)
##     for (pt, value, dist) in ret:
##       echo fmt"point={pt}, value={value}, dist={dist}"
## 
##   # Perform withinRadius searches; notice that the radius distance should use 
##   # the same distance metric as distFunc. The default is squared distance.
##   var ret2 = tree.withinRadius([point.x, point.y], radius=5.0, sortResults=true)
##   for (pt, value, dist) in ret2:
##     echo fmt"point={pt}, value={value}, dist={dist}"
## 
##   # Perform withinRange searches
##   var 
##     min: array[2, float] = [0.0, 0.0]
##     max: array[2, float] = [10.0, 10.0]
##     hyperRect = newHyperRectangle(min, max)
##   
##   var ret = tree.withinRange(hyperRect)
##   for (pt, value) in ret:
##     echo fmt"point={pt}, value={value}"

import algorithm, math

type KdPoint*[K: static[int]] =
      array[K, float]
      ## A KdPoint is a location in K-dimensional space.

# sqrDist returns the square distance between two points.
func sqrDist[K: static[int]](self, other: KdPoint[K]): float =
    result = 0.0
    for i in 0..<K:
        result += (self[i] - other[i]) * (self[i] - other[i])

type DistFunc*[K: static[int]] = 
    proc (x, y: KdPoint[K]): float {.closure.}
    ## The signature of a distance function used to calculate the distance between two KdPoints, 
    ## returning a float. The default DistFunc is `sqrDist`, which returns the squared-distance.

type KdNode[K: static[int]; T] = ref object
    left, right: KdNode[K, T]
    point: KdPoint[K]
    data: T
    splitDimension: int

func newNode[K: static[int]; T](point: KdPoint[K], data: T): KdNode[K, T] =
    new(result)
    result.point = point
    result.data = data

type KdTree*[K: static[int]; T] = object
    ## A k-d tree data structure that allows efficient spatial querying on point distributions.
    ## The tree can hold any generic type `T` data associated with each `K` dimensional point.
    root*: KdNode[K, T]
    len: int
    distFunc: DistFunc[K]

proc buildTree[K: static[int]; T](nodes: var seq[KdNode[K, T]], depth = 0): KdNode[K, T] =
    let numPoints = len(nodes)
    if numPoints > 1:
        let split = depth mod K
        proc kdNodeCmp(x, y: KdNode[K, T]): int =
            if x.point[split] < y.point[split]: -1
            elif x.point[split] == y.point[split]: 0
            else: 1

        nodes.sort(kdNodeCmp)
        let m = (numPoints / 2).int
        result = nodes[m]
        # echo fmt"point={result.point}, data={result.data}"

        result.splitDimension = split
        var left = nodes[0..m-1]
        result.left = buildTree(left, depth+1)
        var right = nodes[m+1..high(nodes)]
        result.right = buildTree(right, depth+1)
    elif numPoints == 1:
        result = nodes[0]
        # echo fmt"leaf point={result.point}, data={result.data}"
    else:
        result = nil

proc newKdTree*[K: static[int]; T](pointData: openArray[(KdPoint[K], T)], distFunc: DistFunc[K] = sqrDist): KdTree[K, T] =
    ## Constructs a k-d tree by bulk-loading an array of point-data tuples, where the associated data is 
    ## of any generic type `T`. Notice that this way of constructing a KdTree should be preferred over 
    ## adding points individually because the resulting tree will be balanced, which will make for more 
    ## efficient search operations. The default `distFunc` is the squared distance, which is returned from
    ## most search functions (except `withinRange`).
    ## 
    ## .. code-block:: nim
    ##  let pointsAndValues = [([2.0, 3.0], 1), 
    ##                        ([5.0, s4.0], 2), 
    ##                        ([9.0, 6.0], 3), 
    ##                        ([4.0, 7.0], 4), 
    ##                        ([8.0, 1.0], 5), 
    ##                        ([7.0, 2.0], 6)]
    ## 
    ##  # Create a 2-D tree of int's
    ##  var tree = newKdTree[2, int](pointsAndValues)
    ##  
    ##  # A custom distance function; the default is a squared distance
    ##  proc myDistFunc(p1, p2: array[2, float]): float {.closure.} =
    ##    result = 0.0
    ##    for i in 0..<len(p1):
    ##      result += (p1[i] - p2[i]) * (p1[i] - p2[i])
    ##    result = sqrt(result)
    ##
    ##  tree = newKdTree[2, int](pointsAndValues, distFunc=myDistFunc)

    doAssert len(pointData) > 0, "Point data appears to be empty"

    var nodes = newSeqOfCap[KdNode[K, T]](len(pointData))
    for p in pointData:
        nodes.add(newNode[K, T](p[0], p[1]))

    result.root = buildTree(nodes)
    result.len = len(nodes)
    result.distFunc = distFunc

proc newKdTree*[K: static[int]; T](points: openArray[KdPoint[K]], data: openArray[T], distFunc: DistFunc[K] = sqrDist): KdTree[K, T] =
    ## Constructs a k-d tree by bulk-loading arrays of points and associated data values of any generic 
    ## type `T`. Notice that this way of constructing a KdTree should be preferred over adding points 
    ## individually because the resulting tree will be balanced, which will make for more efficient 
    ## search operations. The default `distFunc` is the squared distance, which is returned from
    ## most search functions (except `withinRange`).
    ## 
    ## .. code-block:: nim
    ##  let points = [[2.0, 3.0], [5.0, 4.0], [9.0, 6.0], [4.0, 7.0], [8.0, 1.0], [7.0, 2.0]]
    ##  let values = [1, 2, 3, 4, 5, 6]
    ## 
    ##  # Create a 2-D tree of int's
    ##  var tree = newKdTree[2, int](points, values)
    ## 
    ##  # A custom distance function; the default is a squared distance
    ##  proc myDistFunc(p1, p2: array[2, float]): float {.closure.} =
    ##    result = 0.0
    ##    for i in 0..<len(p1):
    ##      result += (p1[i] - p2[i]) * (p1[i] - p2[i])
    ##    result = sqrt(result)
    ##
    ##  tree = newKdTree[2, int](points, values, distFunc=myDistFunc)

    doAssert len(points) == len(data), "Points and data arrays must be the same size."
    doAssert len(points) > 0, "Point data appears to be empty"

    var nodes = newSeqOfCap[KdNode[K, T]](len(points))
    for i in 0..<len(points):
        nodes.add(newNode[K, T](points[i], data[i]))

    result.root = buildTree(nodes)
    result.len = len(nodes)
    result.distFunc = distFunc

func add*[K: static[int]; T](tree: var KdTree[K, T], point: KdPoint[K], data: T) = 
    ## This function can be used to add single points, and their associated data of any 
    ## generic type `T` to an existing KdTree object. Notice that this can result in an 
    ## unbalanced tree which is suboptimal for search operation efficiency. Use the 
    ## `rebalance` function after adding a batch of individual points to the tree.

    var node = newNode(point, data)
    var it = tree.root
    var depth = 0
    while it != nil:
        if node.point[it.splitDimension] <= it.point[it.splitDimension]:
            if it.left == nil:
                node.splitDimension = (depth + 1) mod K
                it.left = node
                return
            it = it.left
        else:
            if it.right == nil:
                node.splitDimension = (depth + 1) mod K
                it.right = node
                return
            it = it.right

        depth += 1

    tree.len += 1

func len*[K: static[int]; T](tree: KdTree[K, T]): int = 
  ## Returns the number of nodes contained within the KdTree.
  tree.len

func height[K: static[int]; T](node: var KdNode[K, T]): int =
    if node == nil:
        return 0

    var lht = node.left.height()
    var rht = node.right.height()
    result = max(lht, rht) + 1
    
func height*[K: static[int]; T](tree: var KdTree[K, T]): int =
    ## Returns the height of the KdTree.
    
    result = height(tree.root)

func isBalanced*[K: static[int]; T](tree: var KdTree[K, T]): int =
    ## Returns the value of the left tree height - right tree height. The larger the 
    ## value magnitude, the more unbalanced the tree is (some say an unbalanced tree 
    ## is any with an absolute magnitude greater than 1). The sign indicates the direction
    ## of skew, with negative values indicating a left-skewed tree and positive values
    ## indicated a right-skewed tree.

    result = height(tree.root.left) - height(tree.root.right)

func rebalance*[K: static[int]; T](tree: var KdTree[K, T]) =
    ## Re-balances an unbalanced KdTree. Note that the original tree structure can be 
    ## completely modified by this function. Use this function after adding a batch
    ## of individual nodes to the tree using the `add` function.

    # collect all the tree's nodes
    var nodes = newSeqOfCap[KdNode[K, T]](len(tree))

    var stack: seq[KdNode[K, T]] = @[tree.root]
    while stack.len > 0:
        var n = stack.pop()
        if n != nil:
            nodes.add(newNode(n.point, n.data))

            stack.add(n.left)
            stack.add(n.right)

    tree.root = buildTree(nodes)
    tree.len = len(nodes)

proc nearestNeighbour*[K: static[int]; T](tree: var KdTree[K, T], point: KdPoint[K]): (KdPoint[K], T, float) =
    ## Returns the nearest neighbour of an input target point, the data associated with the nearest neighbour, 
    ## and the distance between the target point and the nearest neighbour. Notice that the returned distance
    ## uses the distance metric based on the distFunc parameter when the tree is created. By default, and if
    ## unspecified, this metric is the squared distance.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   let (pt, value, dist) = tree.nearestNeighbour([x, y])
    ##   echo fmt"point={pt}, value={value}, dist={dist}"

    var 
        stack: seq[KdNode[K, T]] = @[tree.root]
        minDist: float = Inf
        dist: float
        diff: float
        split: int
    while stack.len > 0:
        var n = stack.pop()
        # echo fmt"{n.point}, {n.data}"
        dist = tree.distFunc(point, n.point)
        if dist < minDist:
            minDist = dist
            result = (n.point, n.data, minDist)

        split = n.splitDimension
        if point[split] <= n.point[split]:
            if n.left != nil:
                stack.add(n.left)
            
            if n.right != nil:
                diff = point[split] - n.point[split]
                if minDist > diff*diff:
                    stack.add(n.right)
            
        else:
            if n.right != nil:
                stack.add(n.right)
            
            if n.left != nil:
                diff = point[split] - n.point[split]
                if minDist > diff*diff:
                    stack.add(n.left)

proc nearestNeighbours*[K: static[int]; T](tree: var KdTree[K, T], point: KdPoint[K], numNeighbours: int): seq[(KdPoint[K], T, float)] =
    ## Returns a specified number (`numNeighbours`) of nearest neighbours of a target point (`point`). Each return point 
    ## is accompanied by the associated data, and the distance between the target and return points. Notice that the 
    ## returned distance uses the distance metric based on the distFunc parameter when the tree is created. By default, 
    ## and if unspecified, this metric is the squared distance.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   let ret = tree.nearestNeighbours([x, y], numNeighbours=5)
    ##   for (pt, value, dist) in ret:
    ##     echo fmt"point={pt}, value={value}, dist={dist}"

    doAssert numNeighbours > 0, "The parameter `numNeighbours` must be larger than zero."

    if numNeighbours == 1:
        return @[nearestNeighbour(tree, point)]

    var 
        stack: seq[KdNode[K, T]] = @[tree.root]
        minDist: float = Inf
        dist: float
        diff: float
        split: int

    result = newSeqOfCap[(KdPoint[K], T, float)](numNeighbours)

    while stack.len > 0:
        var n = stack.pop()
        dist = tree.distFunc(point, n.point)
        if dist <= minDist or len(result) < numNeighbours:
            if len(result) == 0:
                result.add((n.point, n.data, dist))
            else:
                for a in 0..<numNeighbours:
                    if dist <= result[a][2]:
                        result.insert((n.point, n.data, dist), a)
                        if len(result) > numNeighbours:
                            discard result.pop()
                        break
                    elif a == high(result) and len(result) < numNeighbours:
                        result.add((n.point, n.data, dist))
                        break

            minDist = result[high(result)][2] # it's actually the largest min distance

        split = n.splitDimension
        if point[split] < n.point[split]:
            if n.left != nil:
                stack.add(n.left)
            
            if n.right != nil:
                diff = point[split] - n.point[split]
                if minDist > diff * diff:
                    stack.add(n.right)
            
        else:
            if n.right != nil:
                stack.add(n.right)
            
            if n.left != nil:
                diff = point[split] - n.point[split]
                if minDist > diff * diff:
                    stack.add(n.left)

proc withinRadius*[K: static[int]; T](tree: var KdTree[K, T], point: KdPoint[K], radius: float, sortResults=false): seq[(KdPoint[K], T, float)] =
    ## Returns all of the points contained in the tree that are within a specified radius of a target point. By default, the
    ## returned points are in an arbitrary order, unless the `sortResults` parameter is set to true, in which case the return
    ## points will be sorted from nearest to farthest from the target. Notice that the returned distance uses the distance 
    ## metric based on the distFunc parameter when the tree is created. By default, and if unspecified, this metric is the 
    ## squared distance.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   var ret = tree.withinRadius([x, y], radius=5.0, sortResults=true)
    ##   for (pt, value, dist) in ret:
    ##      echo fmt"point={pt}, value={value}, dist={dist}"

    var 
        stack: seq[KdNode[K, T]] = @[tree.root]
        dist: float
        split: int

    result = newSeq[(KdPoint[K], T, float)]()
    
    if radius <= 0:
        return result

    let sqRadius = radius * radius
    
    while stack.len > 0:
        var n = stack.pop()
        dist = tree.distFunc(point, n.point)
        if dist <= sqRadius:
            result.add((n.point, n.data, dist))

        split = n.splitDimension
        if point[split] <= n.point[split]:
            if n.left != nil:
                stack.add(n.left)
            
            if n.right != nil:
                if radius > abs(point[split] - n.point[split]):
                    stack.add(n.right)
            
        else:
            if n.right != nil:
                stack.add(n.right)
            
            if n.left != nil:
                if radius > abs(point[split] - n.point[split]):
                    stack.add(n.left)

    if len(result) == 0:
        return result

    if sortResults:
        proc kdNodeCmp(x, y: (KdPoint[K], T, float)): int =
                if x[2] < y[2]: -1
                elif x[2] == y[2]: 0
                else: 1

        result.sort(kdNodeCmp)

type HyperRectangle*[K: static[int]] = object
    ## A HyperRectangle is used by the withinRange search function to identify multi-deminsional ranges.
    min*: KdPoint[K]
    max*: KdPoint[K]

func newHyperRectangle*[K: static[int]](min: KdPoint[K], max: KdPoint[K]): HyperRectangle[K] =
    ## Creates a new HyperRectangle
    result.min = min
    result.max = max

func withinRange*[K: static[int]; T](tree: var KdTree[K, T], rectangle: HyperRectangle): seq[(KdPoint[K], T)] =
    ## Returns all of the points contained in the tree that are within a target HyperRectangle. 
    ## 
    ## .. code-block:: nim
    ##   var 
    ##     min: array[2, float] = [0.0, 0.0]
    ##     max: array[2, float] = [100.0, 100.0]
    ##     hyperRect = newHyperRectangle(min, max)
    ##   
    ##   var ret = tree.withinRange(hyperRect)
    ##   for (pt, value) in ret:
    ##     echo fmt"point={pt}, value={value}"

    var 
        stack: seq[KdNode[K, T]] = @[tree.root]
        split: int
        withinRange: bool

    result = newSeq[(KdPoint[K], T)]()
    
    for i in 0..<K:
        if rectangle.max[i] <= rectangle.min[i]:
            # it's an ill-formed HyperRectangle
            return result
    
    while stack.len > 0:
        var n = stack.pop()
        withinRange = true
        for i in 0..<K:
            if n.point[i] < rectangle.min[i] or n.point[i] > rectangle.max[i]:
                withinRange = false
                break
        
        if withinRange:
            result.add((n.point, n.data))

        split = n.splitDimension
        if rectangle.min[split] <= n.point[split]:
            if n.left != nil:
                stack.add(n.left)
        
        if rectangle.max[split] >= n.point[split]:
            if n.right != nil:
                stack.add(n.right)

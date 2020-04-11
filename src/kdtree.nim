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
## 
##   for a in 0..<numPoints:
##     x = rand(100_000) / 1000
##     y = rand(100_000) / 1000
##     points.add([x, y])
##     values.add(a)
## 
##   echo fmt"Building tree of {numPoints} random points..."
##   var tree = newKdTree[int](points, values)
## 
##   # Perform nearestNeighour searches
##   let numSearches = 10_000
##   for a in 0..<numSearches:
##     x = rand(100_000) / 1000
##     y = rand(100_000) / 1000
##     let (pt, values, dist) = tree.nearestNieghbour([x, y])
## 
##   # Perform nearestNeighours searches
##   let n = 10
##   for a in 0..<numSearches:
##     x = rand(100_000) / 1000
##     y = rand(100_000) / 1000
##     let ret = tree.nearestNieghbours([x, y], n)
##     for (pt, value, dist) in ret:
##       echo fmt"point={pt}, value={value}, dist={dist}"
## 
##   # Perform withinRadius searches
##   var ret2 = tree.withinRadius([point.x, point.y], radius=5.0, squaredDist=false, sortResults=true)
##   for (pt, value, dist) in ret2:
##     echo fmt"point={pt}, value={value}, dist={dist}"
## 
##   # Perform withinRange searches
##   var 
##     min: array[2, float] = [0.0, 0.0]
##     max: array[2, float] = [100.0, 100.0]
##     hyperRect = newHyperRectangle(min, max)
##   
##   var ret = tree.withinRange(hyperRect)
##   for (pt, value) in ret:
##     echo fmt"point={pt}, value={value}"

import algorithm, math
# import strformat

const K* = 2
  ## K is the dimensionality of the points in this package's K-D trees.

type KdPoint* = 
      array[K, float]
      ## A KdPoint is a location in K-dimensional space.

# SqDist returns the square distance between two points.
func sqDist(self, other: KdPoint): float =
    result = 0.0
    for i in 0..<K:
        result += (self[i] - other[i]) * (self[i] - other[i])

type KdNode[T] = ref object
    left, right: KdNode[T]
    point: KdPoint
    data: T
    splitDimension: int

func newNode[T](point: KdPoint, data: T): KdNode[T] =
    new(result)
    result.point = point
    result.data = data

type KdTree*[T] = ref object
    ## A k-d tree data structure that allows efficient spatial querying on point distributions. The
    ## Current implementation is designed for 2-D point data, although other dimensionality is possible
    ## simply by modifying the const `K`.
    root*: KdNode[T]
    len: int

func buildTree[T](nodes: var seq[KdNode[T]], depth = 0): KdNode[T] =
    let numPoints = len(nodes)
    if numPoints > 1:
        let split = depth mod K
        proc kdNodeCmp(x, y: KdNode[T]): int =
            if x.point[split] < y.point[split]: -1
            elif x.point[split] == y.point[split]: 0
            else: 1

        nodes.sort(kdNodeCmp)
        let m = (numPoints / 2).int
        result = nodes[m]

        result.splitDimension = split
        var left = nodes[0..m-1]
        result.left = buildTree(left, depth+1)
        var right = nodes[m+1..high(nodes)]
        result.right = buildTree(right, depth+1)
    elif numPoints == 1:
        result = nodes[0]
    else:
        result = nil

func newKdTree*[T](points: openArray[KdPoint], data: openArray[T]): KdTree[T] =
    ## Constructs a k-d tree with a bulk-loading of points and associated data values of type `T`.
    ## Notice that this should be the preferred way of constructing a KdTree because the resulting
    ## tree will be balanced, which will make for more efficient search operations.

    if len(points) != len(data):
        raise newException(ValueError, "Points and data arrays must be the same size.")

    var nodes = newSeqOfCap[KdNode[T]](len(points))
    for i in 0..<len(points):
        nodes.add(newNode(points[i], data[i]))

    result = new(KdTree[T])
    result.root = buildTree(nodes)
    result.len = len(nodes)

func add*[T](tree: var KdTree[T], point: KdPoint, data: T) = 
    ## This function can be used to add single points, and their associated data of type `T`
    ## to an existing KdTree object. Notice that this can result in an unbalanced tree which
    ## is suboptimal for search operation efficiency. 

    var node = newNode(point, data)

    # insert the node into the tree
    if tree == nil:
        tree.root = node
    else:
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

func len*[T](tree: KdTree[T]): int = 
  ## Returns the number of nodes contained within the KdTree.
  tree.len

func height[T](node: var KdNode[T]): int =
    if node == nil:
        return 0

    var lht = node.left.height()
    var rht = node.right.height()
    result = max(lht, rht) + 1
    
func height*[T](tree: var KdTree[T]): int =
    ## Returns the height of the KdTree.
    if tree == nil:
        return 0

    result = height(tree.root)

func isBalanced*[T](tree: var KdTree[T]): int =
    ## Returns the value of the left tree height - right tree height. The larger the 
    ## value magnitude, the more unbalanced the tree is (some say an unbalanced tree 
    ## is any with an absolute magnitude greater than 1). The sign indicates the direction
    ## of skew, with negative values indicating a left-skewed tree and positive values
    ## indicated a right-skewed tree.

    result = height(tree.root.left) - height(tree.root.right)

func rebalance*[T](tree: var KdTree[T]) =
    ## Re-balances an unbalanced KdTree. Note that the original tree structure can be 
    ## complete modified by this function. Use this function after adding a significant
    ## number of individual nodes to the tree with the `add` function.

    # collect all the tree's nodes
    var nodes = newSeqOfCap[KdNode[T]](len(tree))

    var stack: seq[KdNode[T]] = @[tree.root]
    while stack.len > 0:
        var n = stack.pop()
        if n != nil:
            nodes.add(newNode(n.point, n.data))

            stack.add(n.left)
            stack.add(n.right)

    tree = new(KdTree[T])
    tree.root = buildTree(nodes)
    tree.len = len(nodes)

func nearestNieghbour*[T](tree: var KdTree[T], point: KdPoint, squaredDist=false): (KdPoint, T, float) =
    ## Returns the nearest neighbour of an input target point, the data associated with the nearest neighbour, and the distance
    ## between the target point and the nearest neighbour. Internally, inter-point distances are calculated as squared-distances. 
    ## These are rooted in the return values, which decreases efficency, unless the squaredDist parameter is set to true.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   let (pt, values, dist) = tree.nearestNieghbour([x, y], squaredDist=false)

    var 
        stack: seq[KdNode[T]] = @[tree.root]
        minDist: float = Inf
        dist: float
        diff: float
        split: int
    while stack.len > 0:
        var n = stack.pop()
        dist = point.sqDist(n.point)
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

    if not squaredDist:
        result[2] = sqrt(result[2])

func nearestNieghbours*[T](tree: var KdTree[T], point: KdPoint, numNeighbours: int, squaredDist=false): seq[(KdPoint, T, float)] =
    ## Returns a specified number (`numNeighbours`) of nearest neighbours of a target point (`point`). Each return point 
    ## is accompanied by the associated data, and the distance between the target and return points. Internally, inter-point
    ## distances are calculated as squared-distances. These are rooted in the return values, which decreases efficency, unless
    ## the squaredDist parameter is set to true.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   let ret = tree.nearestNieghbours([x, y], 5, squaredDist=false)
    ##   for (pt, value, dist) in ret:
    ##     echo fmt"point={pt}, value={values}, dist={dist}"

    if numNeighbours <= 0:
        raise newException(ValueError, "The parameter `numNeighbours` must be larger than zero.")

    if numNeighbours == 1:
        return @[nearestNieghbour(tree, point, squaredDist)]

    var 
        stack: seq[KdNode[T]] = @[tree.root]
        minDist: float = Inf
        dist: float
        diff: float
        split: int

    result = newSeqOfCap[(KdPoint, T, float)](numNeighbours)

    while stack.len > 0:
        var n = stack.pop()
        dist = point.sqDist(n.point)
        if dist < minDist or len(result) < numNeighbours:
            if len(result) == 0:
                result.add((n.point, n.data, dist))
            else:
                for a in 0..<numNeighbours:
                    if a == high(result) and len(result) < numNeighbours:
                        result.add((n.point, n.data, dist))
                        break
                    elif dist < result[a][2]:
                        result.insert((n.point, n.data, dist), a)
                        if len(result) > numNeighbours:
                            discard result.pop()
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

    if not squaredDist:
        for a in 0..<numNeighbours:
            result[a][2] = sqrt(result[a][2])

func withinRadius*[T](tree: var KdTree[T], point: KdPoint, radius: float, squaredDist=false, sortResults=false): seq[(KdPoint, T, float)] =
    ## Returns all of the points contained in the tree that are within a specified radius of a target point. By default, the
    ## returned points are in an arbitrary order, unless the `sortResults` parameter is set to true, in which case the return
    ## points will be sorted from nearest to farthest from the target.
    ## 
    ## .. code-block:: nim
    ##   let x = 100.0
    ##   let y = 25.0
    ##   var ret2 = tree.withinRadius([x, y], radius=5.0, squaredDist=false, sortResults=true)
    ##   for (pt, value, dist) in ret2:
    ##      echo fmt"point={pt}, value={value}, dist={dist}"

    var 
        stack: seq[KdNode[T]] = @[tree.root]
        dist: float
        split: int

    result = newSeq[(KdPoint, T, float)]()
    
    if radius <= 0:
        return result

    let sqRadius = radius * radius
    
    while stack.len > 0:
        var n = stack.pop()
        dist = point.sqDist(n.point)
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
        proc kdNodeCmp(x, y: (KdPoint, T, float)): int =
                if x[2] < y[2]: -1
                elif x[2] == y[2]: 0
                else: 1

        result.sort(kdNodeCmp)

    if not squaredDist:
        for a in 0..high(result):
            result[a][2] = sqrt(result[a][2])

type HyperRectangle* = ref object
    ## A HyperRectangle is used by the withinRange search function to identify multi-deminsional ranges.
    min*: KdPoint
    max*: KdPoint

func newHyperRectangle*(min: KdPoint, max: KdPoint): HyperRectangle =
    ## Creates a new HyperRectangle
    new(result)
    result.min = min
    result.max = max

func withinRange*[T](tree: var KdTree[T], rectangle: HyperRectangle): seq[(KdPoint, T)] =
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
        stack: seq[KdNode[T]] = @[tree.root]
        split: int
        withinRange: bool

    result = newSeq[(KdPoint, T)]()
    
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

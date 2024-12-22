import matplotlib.pyplot as plt
from primitives import *

DRAW=True

def naive_vertical_shoot(p : Point, segs):
    '''returns (aseg, apoint, bseg, bpoint) where apoint is the lowest visible point from AND STRICTLY ABOVE p on segment aseg in segs,
    and bpoint is the highest visible point from p AND STRICTLY BELOW on a segment bseg in segs.
    
    NOTE: If is a vertex of a segment in segs, do not return that segment (as in Project 4);
           instead, return one strictly above/below p.'''

    above_point = None
    above_seg = None

    below_point = None
    below_seg = None

    vertical_line = p.vertical_line_thru()

    for s in segs:
        inters = vertical_line.intersect_segment(s)
        if inters is not None:
            if s.contains_point(inters): # just check
                if inters == p: # return a segment strictly above p
                   continue 
                elif inters.is_above(p):
                    if above_point is None or inters.is_below(above_point):
                        above_seg = s
                        above_point = inters
                elif inters.is_below(p):
                    if below_point is None or inters.is_above(below_point):
                        below_seg = s
                        below_point = inters



    # TODO: Implement this method
    #   Extend your implementation from Project 4

    return above_seg, above_point, below_seg, below_point

def all_vd_queries(outer, inners, DRAW=True):

    '''returns a set of tuples that describe the vertical decomposition of THIS polygon (and not the insides of its inner cycles).

    input: outer is an axis-aligned bounding rectangle, represented by a counter-clockwise-ordered list of Points
              inners is a list of inner cycles of the polygon, each represented by a clockwise-ordered list of Points

    output representation:
        a set of tuples of the form (q, visible_seg, visible_point) where, for each vertex q of the inner cycles:
        - visible_seg is an edge of a DIFFERENT cycle (inner or outer) that is visible vertically (above or below) from q
        - visible_point is the point on visible_seg seen from q
        
        NOTE: If the segment visible (above or below) from q belongs to the same inner cycle as q, do NOT include
         the corresponding tuple in the output, as the segment (q,visible_point) lies inside the inner cycle.

    hint: The constraint that visible_seg lies on a different inner cycle can be checked in linear time, but it can be done in only constant time.
          Consider how the two edges incident to a vertex must appear around that vertex if looking, say, upwards, is to look inside the cycle (in
          which case the visible segment above must belong to the same cycle.)
         
    ASSUMPTIONS:
    - outer is a axis-aligned rectangle that contains all vertices in inners
    - no two vertices in inners have the same y-coordinate'''


    vd_tuples = set()
    segs = set()
    inner_verts = set()

    cycles = [outer] + inners  # both outer, inner
    cycle_ids = {} # need to check if on same cycle
    seg_cyc = {}
    cycle_id = 0

    # STORAGE GUYS
    for cycle in cycles:
        cycle_ids[id(cycle)] = cycle_id
        for i in range(len(cycle)):
            segs.add(Segment(cycle[i - 1], cycle[i]))
            seg_cyc[Segment(cycle[i - 1], cycle[i])] = cycle_id
        cycle_id += 1

    # REAL GUYS
    for cycle in cycles:
        cycle_id = cycle_ids[id(cycle)]
        for i in range(len(cycle)):
            q = cycle[i]

            above_seg, above_point, below_seg, below_point = naive_vertical_shoot(q, segs)
            if above_seg is not None and above_point is not None:
                    if seg_cyc[Segment(above_seg.p1, above_seg.p2)] != cycle_id:
                            mid_point = Segment(q, above_point).midpoint()
                            for inner in inners:
                                inside = False
                                p1x, p1y = inner[0].x(), inner[0].y()
                                for i in range(len(inner) + 1):
                                    p2x, p2y = inner[i % (len(inner))].x(), inner[i % (len(inner))].y()
                                    if min(p1y, p2y) < mid_point.y() <= max(p1y, p2y) and mid_point.x() <= max(p1x, p2x):
                                        if p1y != p2y:
                                            xints = (mid_point.y() - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                                        if p1x == p2x or mid_point.x() <= xints:
                                            inside = not inside
                                    p1x, p1y = p2x, p2y
                            if not inside:
                                vd_tuples.add((q, above_seg, above_point))

            if below_seg is not None and below_point is not None:
                    if seg_cyc[Segment(below_seg.p1, below_seg.p2)] != cycle_id:
                            mid_point = Segment(q, below_point).midpoint()
                            for inner in inners:
                                inside = False
                                p1x, p1y = inner[0].x(), inner[0].y()
                                for i in range(len(inner) + 1):
                                    p2x, p2y = inner[i % (len(inner))].x(), inner[i % (len(inner))].y()
                                    if min(p1y, p2y) < mid_point.y() <= max(p1y, p2y) and mid_point.x() <= max(p1x, p2x):
                                        if p1y != p2y:
                                            xints = (mid_point.y() - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                                        if p1x == p2x or mid_point.x() <= xints:
                                            inside = not inside
                                    p1x, p1y = p2x, p2y
                            if not inside:
                                vd_tuples.add((q, below_seg, below_point))
    if DRAW:
        for s in segs:
            s.draw()
        
        for v in inner_verts:
            v.draw()

        plt.show()

    return vd_tuples


def check_seg_inside_inner(segment, inners):
    mid_point = segment.midpoint()
    for inner in inners:
        inside = False
        p1x, p1y = inner[0].x(), inner[0].y()
        for i in range(len(inner) + 1):
            p2x, p2y = inner[i % (len(inner))].x(), inner[i % (len(inner))].y()
            if min(p1y, p2y) < mid_point.y() <= max(p1y, p2y) and mid_point.x() <= max(p1x, p2x):
                if p1y != p2y:
                    xints = (mid_point.y() - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                if p1x == p2x or mid_point.x() <= xints:
                    inside = not inside
            p1x, p1y = p2x, p2y
    return inside


def get_vertical_decomp(outer, inners, DRAW=False):
    '''returns a representation of the vertical decomposition of the given polygon with holes.

    ASSUMPTIONS
    - outer is a axis-aligned rectangle that contains all vertices in inners
    - no two vertices in inners have the same y-coordinate
    '''

    query_answers = all_vd_queries(outer, inners, DRAW=DRAW) # TODO: Complete this; see above

    segs = set()

    for i in range(len(outer)):
        segs.add(Segment(outer[i], outer[i-1]))

    inner_verts = set()

    for inner in inners:
        for i in range(len(inner)):
            inner_verts.add(inner[i])
            segs.add(Segment(inner[i], inner[i-1]))

    pts_on_segs = { seg:{seg.p1, seg.p2} for seg in segs }
    vert_segs = set()

    for q, vis_seg, vis_p in query_answers:
        vert_segs.add(Segment(q,vis_p))
        pts_on_segs[vis_seg].add(vis_p)

    subdiv_segs = set()
    for seg in segs:
        if seg.is_vertical():
            subdiv_segs.add(seg)

        pts = sorted(pts_on_segs[seg])
        for i in range(len(pts)-1):
            subdiv_segs.add(Segment(pts[i], pts[i+1]))

    subdiv_segs = subdiv_segs.union(vert_segs)
    
    adj = {}
    for seg in subdiv_segs:

        if seg.p1 not in adj:
            adj[seg.p1] = []

        if seg.p2 not in adj:
            adj[seg.p2] = []
            
        adj[seg.p1].append(seg.p2)
        adj[seg.p2].append(seg.p1)
        
    for v in adj:
        adj[v] = sorted( adj[v], key=lambda q: -v.angle(q) )

    visited = set()
    trapezoids = []
    for seg in vert_segs:
        
        first, nxt = seg.top, seg.bottom
        cur = first
        if (cur,nxt) not in visited:

            trap = [first]
            while nxt != first:
                trap.append(nxt)
                visited.add((cur,nxt))
                idx = adj[nxt].index(cur)
                nxt_nxt = adj[nxt][(idx+1)%len(adj[nxt])]

                cur, nxt = nxt, nxt_nxt
                    
            visited.add((cur,nxt))
            trapezoids.append(tuple(trap))

        first, nxt = seg.bottom, seg.top
        cur = first
        if (cur,nxt) not in visited:

            trap = [first]
            while nxt != first:
                trap.append(nxt)
                visited.add((cur,nxt))
                idx = adj[nxt].index(cur)
                nxt_nxt = adj[nxt][(idx+1)%len(adj[nxt])]

                cur, nxt = nxt, nxt_nxt
                    
            visited.add((cur,nxt))
            trapezoids.append(tuple(trap))

    if DRAW:
        for s in subdiv_segs:
            if s.top in inner_verts and s.is_vertical():
                s.draw(color='red')
            elif s.bottom in inner_verts and s.is_vertical():
                s.draw(color='blue')
            else:
                s.draw(color='gray')
            s.p1.draw()
            s.p2.draw()

        for trap in trapezoids:
            triangle = plt.Polygon([p.p() for p in trap], facecolor=('lightgray',0.5),)
            plt.gca().add_patch(triangle)

        plt.show()

    return trapezoids, vert_segs

def construct_roadmap(trapezoids, vert_segs, DRAW=False):
    '''returns a representation of the vertices and edges of the roadmap for the given vertical decomposition.
    
    input: a list of trapezoids, each represented as a tuple of Point objects in counter-clockwise order,
    and a set vert_segs of newly created vertical segments (not existing edges of the polygon) of the vertical
    decomposition of a polygon.
    
    output: a 3-tuple (verts, edges, trap_to_vert_map), where:
        - verts: a set of created roadmap vertices, where each vertex is one of two types:
            the midpoint of a vertical segment of the vertical decomposition or the "center of mass" of each trapezoid
            (its x- and y-coords are the averages of the respective coordinates of the trapezoid's vertices)
        - edges: a set of created roadmap edges (Segment objects) between roadmap vertices of different types;
            that is, every edge connects the center of a trapezoid to its (up to four) vertical edges
        - trap_to_vert_map: a dict whose keys are trapezoids in "trapezoids" and the value is the corresponding vertex in "verts"
    '''

    verts = set()
    edges = set()
    trap_to_vert_map = {}

    trap_verts = set()
    edge_verts = set()
    

    # compute one type of vertices, those on vertical segments
    for seg in vert_segs:
        edge_verts.add(seg.midpoint())

    # compute a roadmap vertiex for trapezoid
    for trap in trapezoids:
        total_x = sum(p.x() for p in trap)
        total_y = sum(p.y() for p in trap)
        center = Point(total_x, total_y, len(trap))
        trap_verts.add(center)
        trap_to_vert_map[trap] = center

        for i in range(len(trap)):
            seg = Segment(trap[i-1], trap[i])
            if seg in vert_segs:
                edges.add(Segment(center, seg.midpoint()))
    
    if DRAW:
        for trap in trapezoids:
            for i in range(len(trap)):
                Segment(trap[i-1],trap[i]).draw()

            triangle = plt.Polygon([p.p() for p in trap], facecolor=('lightgray',0.5),)
            plt.gca().add_patch(triangle)

        for e in edges:
            e.draw(color='red')

        for v in trap_verts:
            v.draw(color='red')
        
        for v in edge_verts:
            v.draw(color='blue')
            

        plt.show()

    verts = trap_verts.union(edge_verts)

    return verts, edges, trap_to_vert_map

def get_containing_trapezoid(trapezoids, q):
    '''returns the trapezoid in trapezoids, represented as a tuple of ccw-oriented points,
    that contains the given point q (or otherwise throws an error)
    
    ASSUMPTION: q lies in the interior of (exactly) one trapezoid'''

    for trap in trapezoids:
        contained = True

        for i in range(len(trap)):
            cur = trap[i]
            prv = trap[i-1]

            if not ccw(prv, cur, q):
                contained = False
                break

        if contained:
            return trap

    return None

def get_path(verts, edges, trap_to_vert_map, s, t, DRAW=False):
    '''returns a sequence of segments in the polygon from s to t, obtained as
    an edge to the center of the trapezoid containing s, then roadmap edges to
    the center of the trapezoid containing t, then to t.

    REQUIREMENT: `networkx` package: "pip install networkx"'''

    trapezoids = trap_to_vert_map.keys()
    import networkx as nx
    G = nx.Graph()

    vmap = {}
    for i,v in enumerate(verts):
        vmap[v] = i
        vmap[i] = v
        G.add_node(i,pos=v.p())

    for e in edges:
        u = vmap[e.p1]
        v = vmap[e.p2]
        G.add_edge(u,v)

    strap = get_containing_trapezoid(trapezoids, s)
    ttrap = get_containing_trapezoid(trapezoids, t)

    if strap is None:
        raise ValueError('Given point s={} not in workspace!'.format(s))
    
    if ttrap is None:
        raise ValueError('Given point t={} not in workspace!'.format(t))
    
    sv = trap_to_vert_map[strap]
    tv = trap_to_vert_map[ttrap]

    path = nx.shortest_path(G, vmap[sv], vmap[tv])
    path_edges = [Segment(vmap[ui],vmap[vi]) for ui,vi in zip(path,path[1:])]

    if DRAW:
        for trap in trapezoids:
            for i in range(len(trap)):
                Segment(trap[i-1],trap[i]).draw()

            triangle = plt.Polygon([p.p() for p in trap], facecolor=('lightgray',0.5),)
            plt.gca().add_patch(triangle)

        for e in edges:
            e.draw(color='red')

        for v in verts:
            v.draw(color='red')

        for e in path_edges:
            e.draw(color='blue')
            e.p1.draw(color='blue')
            e.p2.draw(color='blue')

        Segment(s,sv).draw(color='blue')
        Segment(t,tv).draw(color='blue')
        s.draw(color='lightblue', text='s')
        sv.draw(color='blue')
        t.draw(color='lightblue', text='t')
        tv.draw(color='blue')

        plt.show()

    return path_edges

if __name__=='__main__':

    n = 10

    outer = [
        Point(0,0),
        Point(n,0),
        Point(n,n),
        Point(0,n),
    ]

    # BELOW WAS USED FOR README FIGURES
    # inners = [
    #     [
    #         Point(1,1),
    #         Point(2,2),
    #         Point(5,3),
    #         Point(6,2),
    #     ],
    #     [
    #         Point(3,5),
    #         Point(13,16,2),
    #         Point(9,7),
    #         Point(8,3)
    #     ]
    # ]

    # USE BELOW FOR ANALYSIS
    inners = [
        [
            Point(1,1),
            Point(3,3),
            Point(5,1),
        ],
        [
            Point(2,5),
            Point(3,8),
            Point(4,7),
            Point(6,3)
        ],
        [
            Point(5,9,1),
            Point(8,7,1),
            Point(9,3,1),
            Point(7,5,1),
        ],
    ]

    s = Point(1,12,2)
    t = Point(14,19,2)

    trapezoids, vert_segs = get_vertical_decomp(outer, inners, DRAW=True)
    verts, edges, trap_to_vert_map = construct_roadmap(trapezoids, vert_segs, DRAW=True)
    path = get_path(verts, edges, trap_to_vert_map, s, t, DRAW=True)
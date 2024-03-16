import math
import heapq
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

# Point class to store points and also direction vectors

class Point:
    def __init__(self, x1, y1):
        self.x = x1
        self.y = y1



# Event class to store both site and circle events
# isSiteEvent is a boolean flag that is true if it is a site event and false otherwise
# point stores the point at which event occurs : point itself for site event and bottom of the event circle for circle events
# arc stores the corresponding parabolic arc

class Event:
    def __init__(self, p, isSiteEvent):
        self.point = p
        self.isSiteEvent = isSiteEvent
        self.y = p.y
        self.arc = None
    
    # defined so that we can put it in a max heap
    def __lt__(self, other):
        return self.y > other.y

    def __eq__(self, other):
        return self.y == other.y
    
    # defined so that we can also store these in sets
    def __hash__(self):
        return hash((self.point.x, self.point.y,self.isSiteEvent))



# Class to store Voronoi Edge
# siteOnLeft and Right are the respective voronoi sites on the left and the right of the edge
# start and end obviously store the start and end point of the segments
# neighbour is for two part edges, we connect these at the end of the program
# slope and intercept are calculated using the fact that the edge must lie on the perpendicular bisector of the line joining its corresponding sites

class VoronoiEdge:
    def __init__(self, s, l, r):
        self.startPoint = s
        self.siteOnLeft = l
        self.siteOnRight = r
        self.neighbour = None
        self.endPoint = None
        self.slope = 0
        if r.y - l.y != 0:
            self.slope = -1 * ((r.x - l.x) / (r.y - l.y))
        else:
            self.slope = math.inf #voronoi edge is perpendicular to the line joining its adjacent sites
        self.intercept = s.y - self.slope * s.x


# BLElement Class is a Beachline element
# It stores parabolic arcs if it is a leaf element and edges if it is an internal node
# site is a pointer to parabola focus, edge is a pointer to edge based on the type of node
# circleEvent is a pointer to circle Event 
# parent is a pointer to parent node and left and right are pointers to left and right children respectively 
        
class BLElement: #beachline element class which is either a parabola or an edge
    def __init__(self, l=True, s=None, e=None, c=None, p=None):
        self.is_leaf = l     # boolean value which is true if the node is a leaf node i.e. is a parabola
        self.site = s        # if parabola then pointer to its focus
        self.edge = e        # if edge then pointer to it
        self.circleEvent = c # pointer to circle event (disappearance of arc)
        self.parent = p      # pointer to parent node 

        self._left = None 
        self._right = None

    # sets left child with BLElement p
    def set_left(self, p):
        self._left = p
        p.parent = self

    # sets right child with BLElement p
    def set_right(self, p):
        self._right = p
        p.parent = self

    # for a parabola BL element it gets the immediate left neighbouring parabola for p
    @staticmethod
    def get_left_parabola(p):
        return BLElement.get_closest_left_child(BLElement.get_closest_left_ancestor(p))

    # for a parabola BL element it gets the immediate right neighbouring parabola for p
    @staticmethod
    def get_right_parabola(p):
        return BLElement.get_closest_right_child(BLElement.get_closest_right_ancestor(p))

    # gets the first parent to the left on going right continuously (inverse of closest right child)
    @staticmethod
    def get_closest_left_ancestor(p): 
        parent = p.parent
        current = p
        while current == parent._left :
            if not parent.parent : return 0
            current = parent
            parent = parent.parent
        return parent

    # gets the first parent to the right on going left continuously (inverse of closest left child)
    @staticmethod
    def get_closest_right_ancestor(p): 
        parent = p.parent
        current = p
        while current == parent._right:
            if not parent.parent : return 0
            current = parent
            parent = parent.parent
        return parent

    # gets the rightmost child in the left subtree that is a leaf
    @staticmethod
    def get_closest_left_child(p): 
        if not p:
            return None
        curr = p._left
        while not curr.is_leaf:
            curr = curr._right
        return curr

    # gets the left most child in the right subtree that is a leaf. 
    @staticmethod
    def get_closest_right_child(p): 
        if not p:
            return None
        curr = p._right
        while not curr.is_leaf:
            curr = curr._left
        return curr



# Bounding Box class which stores the width and height of the bounding box of the input points

class BoundingBox:
    def __init__(self, width = 0, height = 0):
        self.width = width
        self.height = height



# Voronoi Class which uses all the previously defined classes to compute the voronoi diagram and store all the necessary information
# sites stores the sites input into the program
# edges stores the edges computed by Fortune's Algorithm
# final edges stores all completed edges
# bounding box stores the bounding box of the input
# beachline stores the root node of the beachline BST 
# sweepline store the y coordinate of the sweep line at any instant of time.
# deleted stores the deleted circle events so we know to ignore them during execution
# points stores the computed points during the course of the algorithm
# event_queue is the event queue storing both circle and site events in a maxheap with priority given on y - coordinates

class Voronoi:
    def __init__(self, sites, bounding_box):
        self.sites = sites
        self.edges = []
        self.final_edges = []
        self.bounding_box = bounding_box
        self.beachline = None
        self.sweepline = None
        self.deleted = set()
        self.points = []
        self.event_queue = []

    # runs Fortune's Algorithm on the input set of Points and returns the set of final edges
    def compute_voronoi(self):
        
        for point in self.sites:
            heapq.heappush(self.event_queue, Event(point,True))

        while self.event_queue:
            event = heapq.heappop(self.event_queue)
            self.sweepline = event.y

            if event in self.deleted: # ensure that the event is not a deleted circle event
                self.deleted.remove(event)
                del event
                continue

            if event.isSiteEvent : # event is a site event, we insert new parabola
                self.insert_parabola(event.point)
            else: # otherwise circle event so we need to remove the arc which got squeezed
                self.remove_parabola(event)

            del event


        self.finish_edges(self.beachline) # clean up dangling edges

        # We go through the set of edges and if the edge has a neighbour we connect the two together 
        for edge in self.edges: 
            if edge.neighbour:
                edge.startPoint = edge.neighbour.endPoint
                if edge.startPoint and edge.endPoint:
                    self.final_edges.append(edge)
                del edge.neighbour
                
        return self.points, self.final_edges

    # insert a parabola given a point p in the event queue
    def insert_parabola(self, p):
        # initially empty beachline so self.beachline is of None type
        if not self.beachline:
            self.beachline = BLElement(s = p, l = True) # initialize the beachline with the first encountered site as root

        # if there is a single parabola in our beach line and we encounter a second point. 
        # if this second point is degenerate with the first point
        if self.beachline.is_leaf and self.beachline.site.y - p.y < 1:
            focus = self.beachline.site
            self.beachline.is_leaf = False
            self.beachline.set_left(BLElement(s = focus))
            self.beachline.set_right(BLElement(s = p))

            start = Point((p.x + focus.x)/2, self.bounding_box.height) #this is also not entirely clear!
            self.points.append(start)

            if(p.x > focus.x): # based on x coordinate we define the left and right sites for the edge
                self.beachline.edge = VoronoiEdge(start, focus, p)
            else:
                self.beachline.edge = VoronoiEdge(start, p, focus)
            
            self.edges.append(self.beachline.edge)
            return
        
        # if more than 1 parabola already in our beachline or 1 parabola and new point is not degenerate

        # we first find the parabola of the beachline directly above the newly insered point by traversing the binary tree
        parabola = self.get_parabola(p.x)
        
        # we check if the parabola has a circle event in which case we add it to the deleted set
        if parabola.circleEvent:
            self.deleted.add(parabola.circleEvent)
            parabola.circleEvent = None

        # get the point of intersection of new parabola with the vertically closest beachline arc
        start = Point(p.x, self.get_y(parabola.site, p.x))
        self.points.append(start)

        #define the left and right edges based on which site is to the left and right respectively 
        el = VoronoiEdge(start, parabola.site, p)
        er = VoronoiEdge(start, p, parabola.site)

        #set the right edge as the neighbour of the left edge and add the left edge to our list of voronoi edges
        el.neighbour = er
        self.edges.append(el)

        #now we convert the intersection node into the right edge
        parabola.edge = er # convert it into an edge, so it is no longer a leaf and the newly formed arcs are the leaves 
        parabola.is_leaf = False
        
        pl = BLElement(s = parabola.site) # left parabola is the left part of the original one 
        pnew = BLElement(s = p) # the new one is the parabola with focus at p on the sweep line
        pr = BLElement(s = parabola.site) # right parabola is the right part of the original one

        parabola.set_right(pr) # right child of the "right edge" is pr 
        parabola.set_left(BLElement(l = False)) # left child we initialize it to be empty edge node and then fill in the contents of the left edge
        parabola._left.edge = el # set it as our el

        parabola._left.set_left(pl) # obviously has the left parabola and the new parabola as the left and right children respectively
        parabola._left.set_right(pnew)

        # now we check for possible circle events in the affected neighbourhoods i.e the pl and pr neighbourhoods
        self.check_circle_event(pl)
        self.check_circle_event(pr)

    # pass the x coordinate of the site, and get the corresponding parabola vertically above the event point
    def get_parabola(self,x0): 
        edge = self.beachline
        x = 0

        # traverse the binary search tree using the x coordinates as keys
        while not edge.is_leaf:
            x = self.get_intersection(edge,self.sweepline) # gets the intersection with inner node (edge)
            if(x > x0):
                edge = edge._left
            else:
                edge = edge._right

        parabola = edge # we have found the leaf node
        return parabola

    # gets the intersection point of the parabolas to the left and right of the edge
    def get_intersection(self, edge, y):

        left = self.beachline.get_closest_left_child(edge) # parabola to the immediate left of the edge - rightmost leaf of left subtree
        right = self.beachline.get_closest_right_child(edge) # parabola to the immediate right of the edge - leftmost leaf of right subtree

        l = left.site # focus of left parabola
        r = right.site # focus of right parabola
         
        # y = 1/2(l.y - y)*(x-l.x)^2 + 1/2(l.y + k) solve for the two parabolae
        d1 = 2.0 * (l.y - y)
        d2 = 2.0 * (r.y - y)

        if l.y == y:
            x_coordinate = l.x
            return x_coordinate
        
        if r.y == y:
            x_coordinate = r.y
            return x_coordinate
        
        a = 1.0 / d1 - 1.0 / d2
        b = (-2.0 *l.x)/ d1 + (2.0 * r.x)/ d2
        c = (y + d1/4 + l.x**2/d1) - (y + d2/4 + r.x**2/d2)

        disc = b * b - 4 * a * c

        #in case the two parabolae have the same y coordinate we exploit symmetry 
        if a == 0:
            x_coordinate = (l.x + r.x)/2
            return x_coordinate
        
        x1 = (-b + math.sqrt(disc)) / (2 * a)
        x2 = (-b - math.sqrt(disc)) / (2 * a)

        if l.y < r.y:
            x_coordinate = max(x1, x2)
        else:
            x_coordinate = min(x1, x2)

        return x_coordinate

    # gets the y coordinate given the x coordinate and the focus p of the parabola
    def get_y(self, p, x):

        d = 2 * (p.y - self.sweepline)
        if d == 0:
            return p.y if p.x == x else self.bounding_box.height
        a = 1 / d
        b = -2 * p.x / d
        c = self.sweepline + d / 4 + p.x**2 / d
        
        return a*x**2 + b*x + c
    
    # removes a parabola when it has a circle event. Pass in the circle event as input.
    def remove_parabola(self, event):
        # we remove a parabola only when it has a circle event. 
        # the associated blelement we know is a parabola
        p = event.arc

        pl = self.beachline.get_left_parabola(p)
        pr = self.beachline.get_right_parabola(p)

        if pl.circleEvent:
            self.deleted.add(pl.circleEvent)
            pl.circleEvent = None

        if pr.circleEvent:
            self.deleted.add(pr.circleEvent)
            pr.circleEvent = None

        squeeze_point = Point(event.point.x, self.get_y(p.site, event.point.x))
        self.points.append(squeeze_point) # add it into the set of voronoi points of the diagram

        # we get the corners of the squeezed parabola and extend them to get the squeeze point
        closest_left_ancestor = self.beachline.get_closest_left_ancestor(p)
        closest_right_ancestor = self.beachline.get_closest_right_ancestor(p)

        closest_left_ancestor.edge.endPoint = squeeze_point
        if closest_left_ancestor.edge.startPoint and closest_left_ancestor.edge.endPoint:
            self.final_edges.append(closest_left_ancestor.edge)

        closest_right_ancestor.edge.endPoint = squeeze_point
        if closest_right_ancestor.edge.startPoint and closest_right_ancestor.edge.endPoint:
            self.final_edges.append(closest_right_ancestor.edge)


        # creates a new edge at the higher predecessor
        ancestor = None
        current = p
        while current != self.beachline:
            current = current.parent
            if current == closest_left_ancestor:
                ancestor = closest_left_ancestor
            if current == closest_right_ancestor:
                ancestor = closest_right_ancestor

        ancestor.edge = VoronoiEdge(squeeze_point, pl.site, pr.site)
        self.edges.append(ancestor.edge)

        # deletes the parabola from the tree and the corresponding corner is also overwritten
        grand_parent = p.parent.parent
        if p.parent._left == p:
            if grand_parent._left == p.parent:
                grand_parent.set_left(p.parent._right)
            if grand_parent._right == p.parent:
                grand_parent.set_right(p.parent._right)
        else:
            if grand_parent._left == p.parent:
                grand_parent.set_left(p.parent._left)
            if grand_parent._right == p.parent:
                grand_parent.set_right(p.parent._left)

        del p.parent

        self.check_circle_event(pl)
        self.check_circle_event(pr)
    
    # cleans up dangling edges and adds their endpoints
    def finish_edges(self, edge):
        if edge.is_leaf: # this means it is a parabola hence not an edge so we simply delete the parabola to free up space
            del edge
            return
        
        bounding_x = 0.0
        # we want to bound the edge inside the bounding box, add an arbitrary constant to ensure that edge extends far enough
        direction = Point(edge.edge.siteOnRight.y - edge.edge.siteOnLeft.y, edge.edge.siteOnLeft.x - edge.edge.siteOnRight.x)
        if direction.x > 0.0:
            bounding_x = max(self.bounding_box.width, edge.edge.startPoint.x + 0.1*self.bounding_box.width) 
        else:
            bounding_x = min(0.0, edge.edge.startPoint.x - 0.1*self.bounding_box.width)

        # compute the endpoint now based on the bounding x coordinate
        end = Point(bounding_x, bounding_x * edge.edge.slope + edge.edge.intercept)
        # add the end point of the edge. 
        edge.edge.endPoint = end

        if edge.edge.startPoint and edge.edge.endPoint:
            self.final_edges.append(edge.edge)

        self.points.append(end)

        # recursively call the children in order to traverse the entire tree
        self.finish_edges(edge._left)
        self.finish_edges(edge._right)

    # checks for circle events on a parabola and its left and right neighbour parabolae
    def check_circle_event(self, blelement):

        # we check circle event on a blelement which is always a parabola

        # gets the left and the right neighbouring parabolae of blelement
        pl = self.beachline.get_left_parabola(blelement) 
        pr = self.beachline.get_right_parabola(blelement)

        # some obvious cases to reject are the ones where pl or pr are null(i.e < 3 parabolae so cant be circle event) 
        # or the two have the same focus which again leads to a case where trivially not a circle event
        if not pl or not pr or pl.site == pr.site:
            return
        
        # call our function to get the intersection of the left border edge and right border edge of the blelement parabola which we passed initially
        centre = self.get_circumcentre(self.beachline.get_closest_left_ancestor(blelement).edge, self.beachline.get_closest_right_ancestor(blelement).edge)

        # if this intersection fails (which is possible if the intersection of the edges doesnt agree with the direction vector of the edges)
        # in that case simply return  
        if not centre:
            return


        # now we compute the Radius of the event circle 
        radius = math.sqrt((pl.site.x - centre.x)**2 + (pl.site.y - centre.y)**2)

        # if the circle is too small to reach the sweep line then this is not a possible circle event and we give up 
        if centre.y - radius >= self.sweepline:
            return
        
        # otherwise if it cuts the sweepline then we definitely have a circle event
        # create a new circle event by setting the siteEvent field as false
        # the site of event is point at which circle event circle touches the directrix 
        circle_event = Event(Point(centre.x, centre.y - radius), False) 

        self.points.append(circle_event.point)

        # we link the arc and the circle event together
        blelement.circleEvent = circle_event
        circle_event.arc = blelement
        
        # add it to the event queue
        heapq.heappush(self.event_queue, circle_event)

    # helper function for check circle event which gets the intersection point of the neighbouring edges of the circle event arc
    def get_circumcentre(self,left, right):

        if left.slope - right.slope == 0:
            return None
        
        # We need to deal with infinite slopes
        if (left.slope == math.inf or left.slope == -math.inf) and (right.slope == math.inf or right.slope == -math.inf):
            return None

        if left.slope == math.inf or left.slope == -math.inf:
            x_intersection = left.startPoint.x
            y_intersection = right.slope * x_intersection + right.intercept

            direction_left = Point(left.siteOnRight.y - left.siteOnLeft.y, left.siteOnLeft.x - left.siteOnRight.x)
            direction_right = Point(right.siteOnRight.y - right.siteOnLeft.y, right.siteOnLeft.x - right.siteOnRight.x)

            if (x_intersection - left.startPoint.x)*direction_left.x < 0 or (y_intersection - left.startPoint.y)*direction_left.y < 0 or (x_intersection - right.startPoint.x)*direction_right.x < 0 or (y_intersection - right.startPoint.y)*direction_right.y < 0:
                return None
            
            centre = Point(x_intersection, y_intersection)
            self.points.append(centre)

            return centre

        if right.slope == math.inf or right.slope == -math.inf:
            x_intersection = right.startPoint.x
            y_intersection = left.slope * x_intersection + left.intercept

            direction_left = Point(left.siteOnRight.y - left.siteOnLeft.y, left.siteOnLeft.x - left.siteOnRight.x)
            direction_right = Point(right.siteOnRight.y - right.siteOnLeft.y, right.siteOnLeft.x - right.siteOnRight.x)

            if (x_intersection - left.startPoint.x)*direction_left.x < 0 or (y_intersection - left.startPoint.y)*direction_left.y < 0 or (x_intersection - right.startPoint.x)*direction_right.x < 0 or (y_intersection - right.startPoint.y)*direction_right.y < 0:
                return None
            
            centre = Point(x_intersection, y_intersection)
            self.points.append(centre)

            return centre
        
        x_intersection = (right.intercept - left.intercept) / ( left.slope - right.slope) # (c2 - c1/ m1 - m2)
        y_intersection = right.slope * x_intersection + right.intercept  # plugging it into the formula for the right edge

        #now we check if this intersection is indeed possible based on the direction of both the lines

        direction_left = Point(left.siteOnRight.y - left.siteOnLeft.y, left.siteOnLeft.x - left.siteOnRight.x)
        direction_right = Point(right.siteOnRight.y - right.siteOnLeft.y, right.siteOnLeft.x - right.siteOnRight.x)

        if (x_intersection - left.startPoint.x)*direction_left.x < 0 or (y_intersection - left.startPoint.y)*direction_left.y < 0 or (x_intersection - right.startPoint.x)*direction_right.x < 0 or (y_intersection - right.startPoint.y)*direction_right.y < 0:
            return None
        
        centre = Point(x_intersection, y_intersection)
        self.points.append(centre)

        return centre
    

if __name__ == "__main__":
    n = int(input())

    sites = []

    max_x = -math.inf
    max_y = -math.inf
    min_x = math.inf
    min_y = math.inf

    for i in range(n):
        x, y = map(float, input().split())
        max_x = max(max_x,x)
        max_y = max(max_y,y)
        min_x = min(min_x,x)
        min_y = min(min_y,y)
        sites.append(Point(x, y))

    bb = BoundingBox(max_x + (max_x - min_x)*0.5, max_y + (max_y - min_y)*0.5)

    voronoi = Voronoi(sites, bb)
    points, edges = voronoi.compute_voronoi()

    # Handling Query Logic
    qx, qy = map(float, input(f"\nEnter positive query point (x,y): ").split())

    # first we find the closest vertex
    min_distance = math.inf
    closest_site = sites[0]
    for site in sites:
        distance = math.sqrt((site.x - qx)**2 + (site.y - qy)**2)
        if distance < min_distance:
            min_distance = distance
            closest_site = site

    print(f"\nThe closest site is ({closest_site.x},{closest_site.y})")
    print("\nThe bounding edges are : ")

    i=1
    printed = set()
        
    for edge in edges:
        #prints if it bounds the query point
        if (edge.siteOnLeft == closest_site or edge.siteOnRight == closest_site) and (edge not in printed):
            if edge.endPoint:
                print(f"Edge {i} : start = ({edge.startPoint.x},{edge.startPoint.y}) , end = ({edge.endPoint.x},{edge.endPoint.y})")
            else:
                print(f"Edge {i} : start = ({edge.startPoint.x},{edge.startPoint.y}) , end = inf")
            printed.add(edge)
            i+=1

    points_new = []
    for point in points:
        points_new.append((point.x,point.y))

    sites_new = []
    for point in sites:
        sites_new.append((point.x,point.y))
        
    fig, ax = plt.subplots()
    for edge in edges:
        #this block of code basically throws out the edge if it passes through a point
        flag = 1
        for point in sites:
            # check if point lies on the edge
            if point.x - edge.startPoint.x != 0:
                m = (point.y - edge.startPoint.y)/(point.x - edge.startPoint.x)
            else:
                m = math.inf
            if(edge.slope == m or (edge.slope == -m and abs(m) == math.inf)):
                flag = 0
                break
        if(flag == 0):
            continue

        if edge.endPoint is not None:
            ax.add_patch(Polygon([[edge.startPoint.x, edge.startPoint.y], [edge.endPoint.x, edge.endPoint.y]], closed=None, fill=None))
        else:
            direction = Point(edge.siteOnRight.y - edge.siteOnLeft.y, edge.siteOnLeft.x - edge.siteOnRight.x)
            x = -direction.y + edge.startPoint.x
            y = direction.x + edge.startPoint.y
            ax.plot([edge.startPoint.x, x], [edge.startPoint.y, y], 'k-')

    x, y = zip(*sites_new)
    ax.plot(x, y, 'ro')

    ax.plot(qx,qy,'bo')

    ax.autoscale()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.show()
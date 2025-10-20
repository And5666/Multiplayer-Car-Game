import pygame, math, time, noise, random, socket, threading, sys, queue, signal, os
from pygame.math import Vector2
import numpy as np

#path to settings
pathToSettings = "settings.txt" 

#open the file for the settings folder given the path 
f = open(os.path.abspath(pathToSettings), "rt")
path = "images/"
lines = []

for line in f:
    lines.append(line.strip())

num = ""
for i in lines[3]:
    if i == ";":
        break
    num += i

f.close()

#initilise pygame and set the clock
pygame.init()

#Game Settings
screenWidth = 1600
screenHeight = 900

fps = int(num)
dt = 1 / fps

clock = pygame.time.Clock() 
pygame.display.set_caption("Multiplayer Car Game")
screen = pygame.display.set_mode((screenWidth, screenHeight), vsync = 1)

class physicsObject:
    def __init__(self, x, y, rotation, width, height, mass, colour):
        #Geometry Variables
        self.position = Vector2(x, y) #starts at x and y poitn given
        self.centreOfMassOffset = Vector2(20,10)
        self.width = width
        self.height = height
        #The list of verticies that make up the shape 
        self.vertices = [
            Vector2(self.position.x - (self.width / 2), self.position.y - (self.height / 2)),
            Vector2(self.position.x + (self.width / 2), self.position.y - (self.height / 2)), 
            Vector2(self.position.x + (self.width / 2), self.position.y + (self.height / 2)), 
            Vector2(self.position.x - (self.width / 2), self.position.y + (self.height / 2))] 
        
        #Linear Componants
        self.mass = mass
        self.velocity = Vector2(0,0)
        self.netForce = Vector2(0,0)
        self.acceleration = Vector2(0,0)
        self.gravity = 0.4
        self.restitution = 0.05

        #Anular Componants  
        self.inertia = (1/12) * self.mass * (self.width ** 2 + self.height **2)
        self.torque = 0  
        self.angularVelocity = 0
        self.angle = 0
        self.rotation = rotation
        self.previousRot = 0

        #set mass and intertia 0 for static objects
        if self.mass == 0:
            self.invMass = 0
            self.invInertia = 0
        else:
            self.invMass = 1/self.mass
            self.invInertia = 1/ self.inertia

        self.colour = colour

    def updatePhysics(self):
        #self.addForce(Vector2(0, self.mass * self.gravity))

        #calculate velocity with force and mass.
        self.acceleration = self.netForce * self.invMass * 3779
        self.velocity += self.acceleration * dt 
        #calculate angular velocity with torque and inertia
        self.angularVelocity += self.torque * 3779 * self.invInertia * dt

        #increase the rotation 
        self.rotation += self.angularVelocity * dt

        #calculate the angle of the shape
        self.angle = (self.previousRot - self.rotation)

        for vertex in self.vertices: #loop through each vertex and rotate + move it
            vertex.x = self.rotateVertex(vertex).x #rotate each vertex
            vertex.y = self.rotateVertex(vertex).y 
                
            vertex += self.velocity * dt  #move each vertex individually 

        #Position calculated with velocity
        self.position += self.velocity * dt   #the position variable is used as the centre of mass. Used as the centre vertex
        
        #set net force and torque back to 0 
        self.netForce = Vector2(0,0)
        self.torque = 0

        self.addForce(Vector2(0, self.mass * self.gravity)) #debuging

        #set current rotation to previous rotation
        self.previousRot = self.rotation
        #self.previousPosition = self.position
        
    def renderObject(self):
        #draw the polygon to the screen 
        pygame.draw.polygon(screen, self.colour, self.vertices)

    def rotateVertex(self, vertex): 
        #move vertex to centre
        vertex -= Vector2(self.position.x, self.position.y)

        #rotate the x and y position using formula below
        tempx = vertex.x #vertices x and y values stored in a temp value so that it dosnt spiral into nothing
        tempy = vertex.y  # this prevents the rotation using the new value of x to caluclate it 

        vertex.x = (tempx * math.cos(self.angle/2)) - (tempy * math.sin(self.angle/2));
        vertex.y = (tempx * math.sin(self.angle/2)) + (tempy * math.cos(self.angle/2));

        #move vertex back to its original position
        vertex += Vector2(self.position.x, self.position.y)

        #return the new position of the vertex
        return vertex
    
    def addForce(self, force):
        self.netForce += force 

    def addTorque(self, torque):
        self.torque += torque 

    def addForceAtPosition(self, force, position):
        #calculate vector from position to centre
        ra = position - (self.position + self.centreOfMassOffset)

        #apply the impulses for obejct A and object B
        #self.addForce(force)
        #self.addTorque((np.cross(ra, -force) )

        self.addForce(force * self.invMass)
        self.angularVelocity += (np.cross(ra, -force) * self.invInertia) *dt
     
class GameManager:
    def __init__(self, car):
        #set the car variable
        self.car = car

        #Camera Settings
        self.cameraOffset = Vector2(-500,50)
        self.cameraPosition = Vector2(0,0)
        self.camSpeed = 10

        #set the static list
        self.staticList = []

        #Terrain Variables
        self.terrainWidth = 33  #width and height of the rect
        self.terrainHeight = 50

        self.terrainLength = 4900  #number of points generated
        self.terrainDensity = 15 #Distance between each point - controls smoothness
        self.terrainXScale = 0.0007 #how streached the x axis is
        self.terrainYScale = 400    #how tall it can be
        self.smallerTerrainYScale = 200
        self.terrainPoints = []    #all points on the terrain

        self.terrainImg = pygame.image.load(path + "grass.png").convert_alpha() #terrain image
        self.terrainImg = pygame.transform.scale(self.terrainImg, (self.terrainWidth, self.terrainHeight/2))  

        #Chunk Variables
        self.chunkLength = 70
        self.chunkList = []
        self.dirtChunkList = []
        self.minDistance = 850
        self.chunk = []
        self.dirtObj = []

        self.fuelList = []
        self.fuelPositions = [500, 19200, 36680]

        #Background Variables
        #image path
        f = open(pathToSettings, "rt")
        bg = ""
        for i in f.readlines()[4]:
            if i == ";":
                break
            bg += i 
        f.close()

        if bg == "True":
            self.showBg = True
        else:
            self.showBg = False

        self.backgroundImg = pygame.image.load(path + "sky.png").convert()
        self.backgroundImg = pygame.transform.scale(self.backgroundImg, (screenWidth, screenHeight))  

        self.bgImages = [] #image list
        self.bgPositions = [] #image positions list
        self.bgPositionsOverlap = [] #second image positions list 
        self.bgSpeed = [0.05, 0.1, 0.3, 0.6] #speed for each image

        #load and set background images images
        for i in range(1,5):
            bgImage = pygame.image.load(path + f"bg{i}.png").convert_alpha()
            bgImage = pygame.transform.scale(bgImage, (screenWidth+500, screenHeight+500)) #scale image slightly larger than screen
            self.bgImages.append(bgImage)

            self.bgPositions.append(Vector2(0,-200))  #append image positions
            self.bgPositionsOverlap.append(Vector2(screenWidth,-200))
        
        #start and end gate images
        self.startPost1 = pygame.image.load(path + "flag-post-1.png").convert_alpha()
        self.startPost2 = pygame.image.load(path + "flag-post-2.png").convert_alpha()

        self.endPost1 = pygame.image.load(path + "flag-post-1.png").convert_alpha()
        self.endPost2 = pygame.image.load(path + "flag-post-2.png").convert_alpha()

    def renderObjects(self):
        pygame.display.set_caption("FPS: " + str(round(clock.get_fps(), 2)))

        #calcualte the offset for the camera
        offset = Vector2(self.cameraPosition.x - screenWidth/2, self.cameraPosition.y - screenHeight/2)

        #render the sky
        screen.blit(self.backgroundImg, Vector2(0,0))
        
        #render the background images
        if self.showBg:
            self.parallaxBackground()

        #render half the start and end gate
        screen.blit(self.startPost1, Vector2(1400,500) - offset)
        screen.blit(self.endPost1, Vector2(70400,500) - offset)

        #select the current chunk
        self.loadChunk()

        #loop throguh each static object
        i = 0
        for obj in self.chunk:
            #set the rotation of the image
            self.tempTerrainImg = pygame.transform.rotate(self.terrainImg, math.degrees(obj.rotation))
            #set the images centre to the middle of it's own self

            #work out how much the images must be offset by
            up = Vector2.normalize(obj.vertices[1] - obj.vertices[2])
            
            terrainRect = self.tempTerrainImg.get_rect(center = (obj.position - offset + obj.width/2 * up))
            #draw grass image to screen
            screen.blit(self.tempTerrainImg, terrainRect)

            #check if its the first or last point in the dirt list
            if i == 0 or i == len(self.dirtObj)-1:
                self.dirtObj[i] = terrainRect.center + self.terrainHeight/4 * -up + screenHeight * 2 * Vector2(0,1) #move point down
            else:
                #set point to the objects point 
                self.dirtObj[i] = terrainRect.center + self.terrainHeight/4 * -up

            #i ++
            i+=1


        #draw the fuel cans
        for can in self.fuelList:
            can.drawCan(offset)

        #draw other players cars
        for player in networkManager.players.values():
            player.renderPlayer(offset)

        #draw the car
        self.car.renderCar(offset)

        #render the other half of the start and end gate
        screen.blit(self.startPost2, Vector2(1420,500) - offset)
        screen.blit(self.endPost2, Vector2(70420,500) - offset)

        #draw the dirt chunk
        pygame.draw.polygon(screen, (61,43,31), self.dirtObj)

        #render the cars Ui
        self.car.renderCarUI()

    def updateObjects(self):
        self.car.updateCar(self.chunk) 
        self.camera()

        self.collideObjects()
        self.car.inputs()   

        #check for collision betwen the car and fuel can 
        for can in self.fuelList:
            if can.collision(self.car.carRect):
                self.fuelList.remove(can)
                self.car.fuelTimer = 0
                self.car.outOfFuel = False

    def loadChunk(self):
        #loop through all chunks
        for i in range(0,len(self.chunkList)):
            #calculate the distance between the car and the middle of the chunk
            distance = abs(self.car.position.x - self.chunkList[i][1].x)
            #check if the distance is low enough
            if distance < self.minDistance:
                if i == 0:
                    #combine the chunks
                    tempList = self.chunkList[i][0] + self.chunkList[i+1][0]
                    tempDirtList = self.dirtChunkList[i] +self.dirtChunkList[i+1]

                    #set the chunks
                    self.chunk = tempList
                    self.dirtObj = tempDirtList
                elif i == len(self.chunkList) -1:
                    #combine the chunks
                    tempList = self.chunkList[i-1][0] + self.chunkList[i][0]
                    tempDirtList = self.dirtChunkList[i-1] + self.dirtChunkList[i] 

                    #set the chunks
                    self.chunk = tempList
                    self.dirtObj = tempDirtList
                else:
                    #combine the chunks
                    tempList = self.chunkList[i-1][0] + self.chunkList[i][0] + self.chunkList[i+1][0]
                    tempDirtList = self.dirtChunkList[i-1] + self.dirtChunkList[i] +self.dirtChunkList[i+1]

                    #set the chunks
                    self.chunk = tempList
                    self.dirtObj = tempDirtList
    
    def parallaxBackground(self):
        #loop through each image
        for i in range(0,4):
            #calculate the positions based on their set speed and the cars velocity 
            self.bgPositions[i].x -= self.bgSpeed[i] * self.car.velocity.x * dt
            self.bgPositions[i].y -= self.bgSpeed[i] * self.car.velocity.y/2 * dt

            self.bgPositionsOverlap[i].x -= self.bgSpeed[i] * self.car.velocity.x * dt
            self.bgPositionsOverlap[i].y -= self.bgSpeed[i] * self.car.velocity.y/2 * dt

            #reset the image if it reaches either side of the screen
            if self.bgPositions[i].x + screenWidth + 500 < 0:
                #if image reaches the left of the screen reset position
                self.bgPositions[i].x = self.bgPositionsOverlap[i].x + screenWidth + 500
            elif self.bgPositions[i].x  > screenWidth + 500:
                #if image reaches the right of the screen reset position
                self.bgPositions[i].x = self.bgPositionsOverlap[i].x - (screenWidth + 500)

            if self.bgPositionsOverlap[i].x + screenWidth + 500 < 0:
                #if image reaches the left of the screen reset position
                self.bgPositionsOverlap[i].x = self.bgPositions[i].x + screenWidth + 500
            elif self.bgPositionsOverlap[i].x  > screenWidth + 500:
                #if image reaches the right of the screen reset position
                self.bgPositionsOverlap[i].x = self.bgPositions[i].x - (screenWidth + 500)

            #draw the images
            tempImg = self.bgImages[i]
            screen.blit(self.bgImages[i], self.bgPositions[i])
            screen.blit(tempImg, self.bgPositionsOverlap[i])
    
    def camera(self):
        #set the camera target
        cameraTarget = self.car.position - self.cameraOffset

        #calculate the distance between the target and its actual position
        camGap = self.cameraPosition - cameraTarget 
        gap = math.sqrt(camGap.x **2 + camGap.y**2)
        #lerp the current position to the target position depending on the gap speed 
        self.cameraPosition.x = self.lerp(self.cameraPosition.x, cameraTarget.x, dt * gap *self.camSpeed)
        self.cameraPosition.y = self.lerp(self.cameraPosition.y, cameraTarget.y, dt * gap *self.camSpeed)

    def spawnTerrain(self):
        #generate points
        #what y level to use of the noise image
        ylevel = 1
        for i in range(0, self.terrainLength):
            #genarate each point
            x = i * self.terrainDensity
            #y level calculated using perlin noise
            smallerY = noise.pnoise1(x  * self.terrainXScale*5, ylevel, repeat = 99999999) * self.terrainYScale/10 + screenHeight/2
            y = noise.pnoise1(x * self.terrainXScale, ylevel, repeat = 99999999) * self.terrainYScale + screenHeight/2 + smallerY

            #set point and append list
            currentPoint = Vector2(x, y)
            self.terrainPoints.append(currentPoint)
        
        #generate physics objects
        for i in range(0,len(self.terrainPoints)):
            #loop through all points in the list and assin current and next point
            currentPoint = self.terrainPoints[i]

            #prevent list going out of bounds
            if i == len(self.terrainPoints) -1:
                continue
            nextPoint = self.terrainPoints[i + 1]

            #work out gradient
            gradient = nextPoint - currentPoint
            #calculate hypotinuse length
            hypotLength = gradient.magnitude()
            #change in height
            oppositeLength = gradient.y

            #sinx = hypot / opposite - trig to work out the angle
            angle =  -math.degrees(math.asin(math.radians(oppositeLength / hypotLength))) 

            #spawn the block at the middle of the two points
            midPoint = Vector2((nextPoint.x + currentPoint.x) /2, (nextPoint.y + currentPoint.y) /2)

            #create the ground block
            groundBlock = physicsObject(midPoint.x , midPoint.y + self.terrainWidth, angle, self.terrainWidth, self.terrainHeight, 0, (255,70,90))

            #append to list and update its rotation
            groundBlock.updatePhysics()

            self.staticList.append(groundBlock)
        
        #split static list into chunks 
        currentChunk = []
        currentDirtChunk = []
        count = 0
        #loop through each item in the static list
        for i in range(0, len(self.staticList)):
            #check if its the start of the chunk
            if count == 0:
                currentDirtChunk.append(self.terrainPoints[i])

            #append the current list
            currentChunk.append(self.staticList[i])
            currentDirtChunk.append(self.terrainPoints[i])

            #check if the max numer of items in a chunk has been reached
            if count == self.chunkLength:
                #define the chunk
                tempChunk = [currentChunk, currentChunk[len(currentChunk) // 2].position]
                currentDirtChunk.append(currentDirtChunk[self.chunkLength-1])

                #add the chunk to the chunk list
                self.chunkList.append(tempChunk)
                self.dirtChunkList.append(currentDirtChunk)

                #reset the current chunk list and counter
                currentDirtChunk = []
                currentChunk = []
                count = 0
        
            #incriment counter
            count+= 1
        
        #spawn in the fuel canisters
        for xPosition in self.fuelPositions:
            fuelPosition = self.raycast(Vector2(xPosition, -10000), Vector2(0,1), 20000, self.staticList)[0]

            currentFuelCan = fuelCan(fuelPosition.x, fuelPosition.y- 50)
            self.fuelList.append(currentFuelCan)
        
    def collideObjects(self):
        #loop through each static object  
        for j in range(0,len(self.chunk)):
            #find the longest side for each object
            maxSideA = max(self.car.width, self.car.height)
            maxSideB = max(self.chunk[j].width, self.chunk[j].height)

            #check if the distance is less than longest sides added
            if self.distanceToObject(self.car, self.chunk[j]) < ((maxSideA + maxSideB) / 2):
                #check for collision if they are close to each other
                self.intersectPolygons(self.car, self.chunk[j])

    def intersectPolygons(self, objectA, objectB): 
        verticesA = objectA.vertices
        verticesB = objectB.vertices

        #set depth to the highest possible value. define depth
        depth = math.inf
        normal = Vector2(0,0)

        for i in range(0,len(verticesA)): #loop through all vertices in shape A
            #set the current vectors 
            vectorA = verticesA[i]  
            vectorB = verticesA[(i + 1) % len(verticesA)] #using mod to make sure I do not exeed array bounds

            #find the gradient of the edge between the two vertexs
            edge = vectorB - vectorA

            #find the normal of the edge
            axis = Vector2(-edge.y, edge.x)
            axis = axis.normalize()

            #set min and max values for A and B using projectVertices() output
            projectVerticesAList = self.projectVertices(verticesA, axis)
            minA = projectVerticesAList[0]
            maxA = projectVerticesAList[1]

            projectVerticesBList = self.projectVertices(verticesB, axis)
            minB = projectVerticesBList[0]
            maxB = projectVerticesBList[1]

            #check if minimum or maximum points have a gap between them
            if minA >= maxB or minB >= maxA:
                #there is a gap so no collision
                return False
            
            #work out the lowest depth 
            axisDepth = min(maxB - minA, maxA - minB)

            #set the lowest depth (if it is the lowest so far) and use the axis from it
            if axisDepth < depth:
                depth = axisDepth
                normal = axis
            
        #same loop but for shape B
        for i in range(0,len(verticesB)): #loop through all vertices in shape B
            #set the current vectors 
            vectorA = verticesB[i]  
            vectorB = verticesB[(i + 1) % len(verticesB)] #using mod to make sure I do not exeed array bounds

            #find the gradient of the edge between the two vertexs
            edge = vectorB - vectorA
            #find the normal of the edge
            axis = Vector2(-edge.y, edge.x)
            axis = axis.normalize()

            #set min and max values for A and B using projectVertices() output
            projectVerticesAList = self.projectVertices(verticesA, axis)
            minA = projectVerticesAList[0]
            maxA = projectVerticesAList[1]

            projectVerticesBList = self.projectVertices(verticesB, axis)
            minB = projectVerticesBList[0]
            maxB = projectVerticesBList[1]

            #check if minimum or maximum points have a gap between them
            if minA >= maxB or minB >= maxA:
                #there is a gap so no collision
                return False
            
            #work out the lowest depth 
            axisDepth = min(maxB - minA, maxA - minB)

            #set the lowest depth (if it is the lowest so far) and use the axis from it
            if axisDepth < depth:
                depth = axisDepth
                normal = axis

        #collision found 
        depth /= normal.length() * dt
        normal = normal.normalize()
        #pygame.draw.line(screen, (0,255,0), objectA.position -100 * normal, objectA.position +100 * normal) #debugging

        direction = objectB.position - objectA.position

        if np.dot(direction, normal) < 0:
            normal = -normal
        
        self.collision(objectA, objectB, normal, depth, self.findContactPoint(objectA.vertices, objectB.vertices))

    def projectVertices(self, vertices, axis):
        #project a list of vertices onto an axis and return min and max
        #set the min and max as the largest possibe values they can be but its flipped
        minimum = math.inf
        maximum = -math.inf

        #loop through each vertex provided
        for i in range(0,len(vertices)):
            #set the current vertex 
            v = vertices[i]
            
            #project v onto the axis
            projection = np.dot(v, axis)

            #compare minimum and maximum to current minimum and maximum
            if projection < minimum:
                minimum = projection #update values
            if projection > maximum:
                maximum = projection

        #pack values into a list 
        minmax = [minimum, maximum]
        #return values
        return minmax

    def findContactPoint(self, verticesA, verticesB):
        #set defult values for outputs
        contact1 = Vector2(0,0)
        contact2 = Vector2(0,0)

        minDistanceSquared = math.inf

        #loop through all vertices of shape A
        for i in range(0, len(verticesA)):
            #set current vertex
            p = verticesA[i]

            #loop through all edges of shape B
            for j in range(0, len(verticesB)):
                #set the two points for the current edge in shape B
                va = verticesB[j]
                vb = verticesB[(j+1) % len(verticesB)]

                #find the cosest point to the vertex 
                closestPointToVertexList = self.closestVertexToPoint(p, va, vb)

                #unpack the output
                distanceSquared = closestPointToVertexList[0]
                contactPoint = closestPointToVertexList[1]

                if self.closeEnough(distanceSquared, minDistanceSquared):
                    #two contact points found
                    if not (self.closeEnough(contactPoint.x, contactPoint.y) and self.closeEnough(contact1.x, contact1.y)):
                        #making sure the contact point is not the same from both shapes
                        contact2 = contactPoint
                        #find the average of the two points
                        contact = Vector2((contact1.x +contact2.x) /2, (contact1.y + contact2.y) /2) 

                if distanceSquared < minDistanceSquared:
                    #only 1 contact Point 
                    minDistanceSquared = distanceSquared
                    contact1 = contactPoint
                    contact = contact1
        
        #loop through all vertices of shape B
        for i in range(0, len(verticesB)):
            #set current vertex
            p = verticesB[i]

            #loop through all edges of shape A
            for j in range(0, len(verticesA)):
                #set the two points for the current edge in shape A
                va = verticesA[j]
                vb = verticesA[(j+1) % len(verticesA)]

                closestPointToVertexList = self.closestVertexToPoint(p, va, vb)

                distanceSquared = closestPointToVertexList[0]
                contactPoint = closestPointToVertexList[1]

                if self.closeEnough(distanceSquared, minDistanceSquared):
                    #two contactPoints found
                    if not (self.closeEnough(contactPoint.x, contactPoint.y) and self.closeEnough(contact1.x, contact1.y)):
                        #making sure the contact point is not the same from both shapes
                        contact2 = contactPoint
                        #find the average of the two points
                        contact = Vector2((contact1.x +contact2.x) /2, (contact1.y + contact2.y) /2) 

                if distanceSquared < minDistanceSquared:
                    #only 1 contact Point 
                    minDistanceSquared = distanceSquared
                    contact1 = contactPoint
                    contact = contact1
        
        pygame.draw.rect(screen, (255,0,0), pygame.Rect(contact.x,contact.y,10,10))

        #return the contact points and the number of contacts
        return contact

    def closestVertexToPoint(self, p, a, b):
        #line
        ab = b - a
        #vector to point from vertex A
        ap = p - a

        #project ap onto ab
        proj = np.dot(ap, ab)
        #calculate the lenght squared and divied the projection by it
        d = proj / (ab.x**2 + ab.y**2)

        #determine the closest point
        if d <= 0:
            closestPoint = a
        elif d >= 1:
            closestPoint = b
        else:
            closestPoint = a + ab * d

        #calculate the length
        l = closestPoint - p
        #length squared
        lSquared = l.x**2 + l.y**2

        #return the length squared and the closest point
        return [lSquared, closestPoint]

    def closeEnough(self, a, b):
        # returns true of the amount is less than 1/2 a mm (close enough)
        return abs(a-b) < 0.0005

    def collision(self, objectA, objectB, normal, depth, contactPoint):
        if objectA.mass == 0 and objectB.mass == 0:
            #both objects are static so no need to resolve the collision
            return
        #resolve collision
        #calculate the resitution
        e = min(objectA.restitution, objectB.restitution)

        #vectors from centre of object to the contact point
        ra = contactPoint - objectA.position
        rb = contactPoint - objectB.position

        #perpendicular vector of ra and rb
        raPerp = Vector2(-ra.y, ra.x)
        rbPerp = Vector2(-rb.y, rb.x)

        #find the velocity at the contact point
        contactPointVelA = raPerp * (objectA.angularVelocity)
        contactPointVelB = rbPerp * (objectB.angularVelocity)

        #calculate the overall velocity of the contact point
        relativeVelocity = (objectB.velocity - contactPointVelB) - (objectA.velocity - contactPointVelA) 

        #make sure the obejcts are moving towards each other
        velAlongNormal = np.dot(relativeVelocity, normal)
        if velAlongNormal > 0:
            return
        
        #project perpendicular contact vectors onto the normal
        raPerpDotN = np.dot(raPerp, normal)
        rbPerpDotN = np.dot(rbPerp, normal)

        # the denominator of the formula
        denominator = (objectA.invMass + objectB.invMass) + ((raPerpDotN * raPerpDotN ) * objectA.invInertia) + ((rbPerpDotN * rbPerpDotN) * objectB.invInertia) 

        #the scalar value of the impulse
        j = (-(1 + e) * velAlongNormal) / denominator / dt

        #vector value of impulse
        impulse = normal * j 

        #apply the impulses for obejct A and object B
        objectA.velocity += -impulse* objectA.invMass * dt
        objectA.angularVelocity += (np.cross(ra, impulse) * objectA.invInertia) * dt
            
        objectB.velocity += impulse * objectB.invMass * dt
        objectB.angularVelocity +=(-np.cross(rb, impulse) * objectB.invInertia) * dt

        #Move objects out of each other
        percent = 0.5 #usually between .2 and .8
        limit = 0.01 # usually 0.01 to 0.1 (helps to reduce stutter)

        #determines if and where the object should be moved 
        correction =(max(depth - limit, 0) / (objectA.invMass + objectB.invMass)) * percent * -normal  
        
        if objectB.mass == 0:
            #only object B is static so only move object A
            objectA.position += objectA.invMass * correction * dt  
            for vertexA in objectA.vertices:
                #move all vertices the same amount as the centre of mass
                vertexA += objectA.invMass * correction * dt
        elif objectA.mass == 0:
            #only object A is static so only move object A
            objectB.position -= objectB.invMass * correction * dt  
            for vertexB in objectB.vertices:
                #move all vertices the same amount as the centre of mass
                vertexB -= objectB.invMass * correction * dt

    def distanceToObject(self, objectA, objectB):
        #calculate the vector from a to b
        ab = objectB.position - objectA.position

        #return the distance 
        return math.sqrt(ab.x **2 + ab.y **2)

    def lerp(self, current, goal, time):
        #calculate the gap 
        difference = goal - current

        #move towards the target 
        if difference > time:
            return current + time
        if difference < -time:
            return current - time
        
        #return the goal if the same
        return goal

    def signalHandler(signum, frame):
        threadingEvent.set()
    
    def raycast(self, startPoint, direction, length, staticList):
        endPoint = startPoint + (direction * length)
        #set default value for max
        maxLength = math.inf

        #set current contact point
        contactPoint = Vector2(0,0)

        for j in range(0, len(staticList)):
            verticesB = staticList[j].vertices

        #check each edge in shape B against the edge used in shape A
            for i in range(0, len(verticesB)):
                    #set the first and second point to make an edge in shape B
                    vertexB1 = verticesB[i]
                    vertexB2 = verticesB[(i +1) % len(verticesB)]

                    #compare lines and return true if they intersect 
                    if self.intersectSegments(startPoint, endPoint, vertexB1, vertexB2):
                        #calculate gradients
                        
                        #make sure the gradients are not the same
                        if endPoint.x - startPoint.x == 0 : #prevent divide by 0 error
                            #set x
                            x = startPoint.x
                            #set gradients
                            mB = (vertexB2.y - vertexB1.y) / (vertexB2.x - vertexB1.x)
                            cB = vertexB1.y - (mB * vertexB1.x)
                            #set y
                            y = (mB * x) + cB
                        elif vertexB2.x - vertexB1.x == 0:
                            #set x
                            x = vertexB1.x 
                            #set gradients
                            mA = (endPoint.y - startPoint.y) / (endPoint.x - startPoint.x)
                            cA = startPoint.y - (mA * startPoint.x)
                            #set y
                            y = (mA * x) + cA
                        else:
                            #calculate gradients
                            mA = (endPoint.y - startPoint.y) / (endPoint.x - startPoint.x)
                            mB = (vertexB2.y - vertexB1.y) / (vertexB2.x - vertexB1.x)
                        
                            #calculate y intercepts
                            cA = startPoint.y - (mA * startPoint.x)
                            cB = vertexB1.y - (mB * vertexB1.x)

                            #find the x and y
                            x = (cB - cA) / (mA - mB)
                            y = (mA * x) + cA

                        #calculate length 
                        ab = startPoint - Vector2(x,y)
                        length = math.sqrt((ab.x **2) + (ab.y **2))
                        
                        #find the point closest to the car
                        if length < maxLength:
                            maxLength = length
                            contactPoint = Vector2(x,y)

                        self.groundNormal = (vertexB2 - vertexB1).normalize()

        if contactPoint != Vector2(0,0):
            output = [contactPoint, True]
            return output

        output = [endPoint, False]
        return output
    
    def intersectSegments(self, a,b,c,d):
        #either acd should be counter-clockwise or bcd but not both. same for abc and abd 
        #returns true if two segments intersect
        return self.ccwCheck(a,c,d) != self.ccwCheck(b,c,d) and self.ccwCheck(a,b,c) != self.ccwCheck(a,b,d)

    def ccwCheck(self, a,b,c):
        #check if a, b and c are in a counter clockwise direction. returns true if they are
        return (c.y-a.y)*(b.x-a.x) > (b.y-a.y)*(c.x-a.x)
    
    def exitGame():
        signal.signal(signal.SIGINT, GameManager.signalHandler)
        pygame.quit()
        sys.exit()

class Car(physicsObject):
    def __init__(self, x, y, rotation, width, height, mass, colour):
        super().__init__(x, y, rotation, width, height, mass, colour)
        #Input Settings
        self.inputHorizontal = 0

        #Direction Vectors
        self.foward = Vector2(0,0)
        self.up = Vector2(0,0)
        self.groundNormal = Vector2(0,0)
        self.direction = 0

        self.centreOfMassOffset = Vector2(5,20)

        #Car Settings
        self.airDrag = 0.1
        self.angularDrag = 100
        self.maxAngularVelocity = 1.3
        self.carTorque = 2000
        self.downForce = 0

        #Fuel Settings
        self.fuelTime = 100 #in seconds
        self.fuelTimer = 0
        self.outOfFuel = False
        self.fuelNeedleOffset = Vector2(0, 30)

        self.fuelGuage = pygame.image.load(path + "fuel-gauge.png").convert_alpha()
        
        self.fuelIcon = pygame.image.load((path + "fuel-icon.png")).convert_alpha()
        self.fuelWarning = pygame.image.load(path + "fuel-warning.png").convert_alpha()
        self.fuelNeedle = pygame.image.load(path + "meter-needle.png").convert_alpha()
        self.fuelNeedle = pygame.transform.scale(self.fuelNeedle, (32/1.5, 128/1.5))

        #Engine Settings
        self.maxRpm = 6500
        self.engineCurveEnd = 2
        self.maxEngineTorque = 100

        self.engineRpm = 0
        self.gearRatio = 7 #lower number will be faster to accelerate but will have a lower top speed

        #Wheel Settings
        self.wheelRadius = 20
        self.wheelInertia = 3000

        self.brakeForce = 6000
        self.rollingResistance = 0.05
        self.grip = 70

        self.awd = True

        #Pacejka Parameters
        self.stiffness = 12  #4-12
        self.shape = 1.9     #1-2
        self.peak = 0.8     #0.1 - 1
        self.curvature = 0.97  #-10 - 1

        #Suspension settings
        self.springTravel = 35
        self.springStiffness = 45000
        self.damperStiffness = 5000

        #Front Wheel
        self.frontAngularVelocity = 0
        self.frontAngularAcceleration = 0
        self.frontVelocity = Vector2(0,0)
        self.frontRotation = 0

        self.frontSuspensionPoint = Vector2(0,0)  
        self.frontContactPoint = Vector2(0,0)
        self.frontSpringLength = self.springTravel

        self.frontLongitudinalForce = 0
        self.weightFront = 0
        self.frontGrounded = False

        #Back Wheel
        self.backAngularVelocity = 0
        self.backAngularAcceleration = 0
        self.backVelocity = Vector2(0,0)
        self.backRotation = 0

        self.backSuspensionPoint = Vector2(0,0)
        self.backContactPoint = Vector2(0,0)
        self.backSpringLength = self.springTravel

        self.backLongitudinalForce = 0
        self.weightBack = 0
        self.backGrounded = False

        #Car Image Variables
        self.carImg = pygame.image.load(path + "Car.png").convert_alpha()
        self.carImg = pygame.transform.scale(self.carImg, (self.width, self.height+20)) 
        self.carRect = self.carImg.get_rect(center = self.position +10 * self.up)

        self.frontWheelImg = pygame.image.load(path + "tire.png").convert_alpha()
        self.frontWheelImg = pygame.transform.scale(self.frontWheelImg, (self.wheelRadius * 2, self.wheelRadius *2))

        self.backWheelImg = pygame.image.load(path + "tire.png").convert_alpha()
        self.backWheelImg = pygame.transform.scale(self.backWheelImg, (self.wheelRadius * 2, self.wheelRadius *2))

        #Pedal Image Variabels
        self.gasPedal = pygame.image.load(path + "pedal-gas-normal.png").convert_alpha()
        self.gasPedalPressed = pygame.image.load(path + "pedal-gas-pressed.png").convert_alpha()

        self.brakePedal = pygame.image.load(path + "pedal-brake-normal.png").convert_alpha()
        self.brakePedalPressed = pygame.image.load(path + "pedal-brake-pressed.png").convert_alpha()

        #Rpm Meter Variables
        self.tachometer = pygame.image.load(path + "tachometer.png").convert_alpha()
        self.rpmNeedle = pygame.image.load(path + "meter-needle.png").convert_alpha()
        self.rpmNeedle = pygame.transform.scale(self.rpmNeedle, (32/1.5, 128/1.5))
        self.rpmNeedleOffset = Vector2(0, 30)
        self.rpmFont = pygame.font.SysFont("Arial", 30, bold = True)

        #Game / ui Settings
        self.startCountdown = False
        self.startGame = False
        self.countdownTimer = 3
        self.waitingForPlayersFont = pygame.font.SysFont("Arial", 80, bold = False)
        self.countdownFont = pygame.font.SysFont("Bahnschrift", 150, bold = True)
        self.timeFont = pygame.font.SysFont("Arial", 40, bold = True)

        self.timer = 0

    def renderCar(self, offset):
        #rotate image
        carImg = pygame.transform.rotate(self.carImg, math.degrees(self.rotation))

        frontWheelImg = pygame.transform.rotate(self.frontWheelImg, math.degrees(self.frontRotation))
        backWheelImg = pygame.transform.rotate(self.backWheelImg, math.degrees(self.backRotation))

        #set the images centre to the middle of it's own self
        self.carRect = carImg.get_rect(center = self.position +10 * self.up - offset )
        screen.blit(carImg, self.carRect)

        #find centre point of wheel
        frontWheelRect = frontWheelImg.get_rect(center = self.frontContactPoint - offset  - self.wheelRadius * -self.up)
        backWheelRect = backWheelImg.get_rect(center = self.backContactPoint - offset - self.wheelRadius* -self.up)

        #draw the wheels
        screen.blit(frontWheelImg, frontWheelRect)
        screen.blit(backWheelImg, backWheelRect)

    def renderCarUI(self):
        self.fuel()
        self.rpmMeter()
        if self.startGame:
            self.timeText()

        #draw the pedals
        #gas pedal
        if self.inputHorizontal == 1:
            screen.blit(self.gasPedalPressed, Vector2(1300, 650))
        else:
            screen.blit(self.gasPedal, Vector2(1300, 650))
        #brake pedal
        if self.inputHorizontal == -1:
            screen.blit(self.brakePedalPressed, Vector2(50, 650))
        else:
            screen.blit(self.brakePedal, Vector2(50, 650))

        if self.startCountdown:
            self.countdown()
        elif not self.startGame:
            text = self.waitingForPlayersFont.render("Waiting For Players...", True, (0,0,0))
            screen.blit(text, (screenWidth/2 -300, screenHeight/2 - 300))

    def updateCar(self, staticList):
        self.updatePhysics()
        #update the bodys physics
        self.weightTransfer()

        self.forward = Vector2.normalize(self.vertices[1] - self.vertices[0])
        self.up = -Vector2(-self.forward.y, self.forward.x)

        self.carPhysics()
        self.wheelPhysics(staticList)

        networkManager.sendPlayerInformation(self.position, self.rotation, self.up, self.frontContactPoint, self.backContactPoint, self.velocity.x)

    def carPhysics(self):
        #calculate rolling velocity
        if self.frontGrounded:
            self.frontAngularVelocity = (self.velocity.magnitude() / self.wheelRadius /2) * self.direction 
        if self.backGrounded:
            self.backAngularVelocity = (self.velocity.magnitude() / self.wheelRadius /2) * self.direction

        #find the wheel rpm from the angular velocity
        frontWheelRpm = abs(self.frontAngularVelocity / (2 * math.pi/60))
        rearWheelRpm = abs(self.backAngularVelocity / (2 * math.pi/60))
        
        #calculate rpm of engine
        self.engineRpm = self.lerp(self.engineRpm, ((rearWheelRpm + frontWheelRpm)/2 * self.gearRatio * 3), 30000 * dt )

        #make the engine idle at 1000 rpm and redline at the maxRpm
        if not self.outOfFuel:
            self.engineRpm = pygame.math.clamp(self.engineRpm, 750, self.maxRpm)
        else:
            self.engineRpm = 0

        #calculate the engine torque given the rpm
        x = self.engineRpm / (self.maxRpm + 500)
        engineTorque = 0
        if x != 0 and not self.outOfFuel:
            engineTorque = self.maxEngineTorque * (-4 * (x-0.5) **2 + 1)
            #engineTorque = self.maxEngineTorque *((1/(x * math.sqrt(2* math.pi))) * math.e **(-0.5*(np.log(x)) **2))
        
        #set the direction the car is moving in
        driveTorque = 0
        frontBrakeTorque = 0
        backBrakeTorque = 0
        if self.startGame:
            if self.velocity.x > 1:
                #moving forwards
                self.direction = 1
                if self.inputHorizontal < 0:
                    #if moving forwards and pressing backwards, apply the brake
                    frontBrakeTorque = self.wheelRadius * self.brakeForce * -self.frontAngularVelocity
                    backBrakeTorque = self.wheelRadius * self.brakeForce * -self.backAngularVelocity 

                else:
                    driveTorque = engineTorque * self.gearRatio * self.inputHorizontal * 1000
            elif self.velocity.x < 1:
                #moving backwards
                self.direction = -1
                driveTorque = engineTorque * self.gearRatio * self.inputHorizontal * 1000
        else:
            if self.velocity.x > -1:
                #moving forwards
                self.direction = 1
            else:
                self.direction = -1
            frontBrakeTorque = self.wheelRadius * self.brakeForce * 10 * -self.frontAngularVelocity
            backBrakeTorque = self.wheelRadius * self.brakeForce * 10 * -self.backAngularVelocity 

        #calculate the rolling torque supplied to the wheels
        if self.frontGrounded:
            frontRollingTorque = self.wheelRadius * self.rollingResistance * self.frontVelocity.magnitude() * -self.frontAngularVelocity
        else:
            frontRollingTorque = self.wheelRadius * self.rollingResistance * 5 * self.frontVelocity.magnitude() * -self.frontAngularVelocity
        if self.backGrounded:
            backRollingTorque = self.wheelRadius * self.rollingResistance * self.backVelocity.magnitude() * -self.backAngularVelocity
        else:
            backRollingTorque = self.wheelRadius * self.rollingResistance * 5 * self.backVelocity.magnitude() * -self.backAngularVelocity

        if self.awd:    
            #calculate the total torque supplied to the wheels
            frontTotalTorque = driveTorque/2 + frontBrakeTorque + frontRollingTorque
            backTotalTorque = driveTorque/2 + backBrakeTorque + backRollingTorque
        else:
            #calculate the total torque supplied to the wheels
            frontTotalTorque = frontBrakeTorque + frontRollingTorque
            backTotalTorque = driveTorque + backBrakeTorque + backRollingTorque

        #apply the torque to the wheels
        self.frontAngularAcceleration = frontTotalTorque / self.wheelInertia
        self.frontAngularVelocity += self.frontAngularAcceleration * dt 

        self.backAngularAcceleration = backTotalTorque / self.wheelInertia
        self.backAngularVelocity += self.backAngularAcceleration * dt 

        #calculate the slip of the wheel
        frictionCoefList = self.slipRatio()

        #calculate longitudinal force 
        self.frontLongitudinalForce = frictionCoefList[0] * self.grip * self.weightFront
        self.backLongitudinalForce = frictionCoefList[1] * self.grip * self.weightBack

        if self.frontGrounded or self.backGrounded:
            #if both grounded add downforce
            #calculate downforce
            downforce = self.downForce * (1 +self.velocity.magnitude())
            #calmp the downforce
            downforce = pygame.math.clamp(downforce, self.downForce, 20000)

            #apply the downforce
            self.addForce(downforce * self.velocity.magnitude() * -self.up  * self.invMass )

            #apply air resistance
            self.addForce(self.airDrag * (Vector2(self.velocity.x **2, self.velocity.x **2) * -self.direction)  * self.invMass )

        if not (self.frontGrounded or self.backGrounded):
            #apply a resistive force if the angular velocity is too high 
            resistiveTorque = self.angularDrag * -self.angularVelocity 
            #calculate spin torque on car
            carTorque = self.carTorque * self.inputHorizontal 

            if -self.angularVelocity > 0 and -carTorque < 0:
                self.addTorque(carTorque + (-self.angularVelocity * 10) /(3779 * self.invInertia))
            elif -self.angularVelocity < 0 and -carTorque > 0:
                self.addTorque(carTorque - (self.angularVelocity * 10) /(3779 * self.invInertia))
            else:
                if abs(self.angularVelocity) < self.maxAngularVelocity:
                    self.addTorque(carTorque)

            self.addTorque(resistiveTorque)
        
        
        
        #print((self.velocity.magnitude()/34))
        #print(self.engineRpm)
        #print(engineTorque)

    def wheelPhysics(self, staticList):
        #set suspsnsion points
        self.frontSuspensionPoint = self.vertices[2] - (28.5 * self.forward ) + 25 * self.up
        self.backSuspensionPoint = self.vertices[3] + (26.5 * self.forward ) + 25 * self.up

        #set spring lengths
        frontSuspensionList = self.suspensionForce(self.frontSuspensionPoint, self.frontSpringLength, staticList)
        self.frontSpringLength = frontSuspensionList[1]
        backSuspensionList = self.suspensionForce(self.backSuspensionPoint, self.backSpringLength, staticList)
        self.backSpringLength = backSuspensionList[1]

        #set contactPoints
        self.frontContactPoint = frontSuspensionList[2]
        self.backContactPoint = backSuspensionList[2]

        #apply suspension force
        self.addForceAtPosition(frontSuspensionList[0], self.frontSuspensionPoint)
        self.addForceAtPosition(backSuspensionList[0], self.backSuspensionPoint)

        #calculate the speed at the wheels
        rf = self.frontContactPoint - self.frontSuspensionPoint
        rb = self.backContactPoint - self.backSuspensionPoint
        contactPointVelF = rf * -self.angularVelocity
        contactPointVelB = rb * -self.angularVelocity
        #set speed at wheels
        self.frontVelocity = self.velocity + contactPointVelF
        self.backVelocity = self.velocity + contactPointVelB

        
        #check if front wheel is grounded
        if frontSuspensionList[3]:
            #set grounded to true
            self.frontGrounded = True
            
            #apply the longitudinal force 
            self.addForceAtPosition(self.frontLongitudinalForce * self.groundNormal, self.frontContactPoint)
        else:   
            self.frontGrounded = False

        #check if back wheel is grounded
        if backSuspensionList[3]:
            #set grounded to true
            self.backGrounded = True
            
            #apply the longitudinal force 
            self.addForceAtPosition(self.backLongitudinalForce * self.groundNormal, self.backContactPoint)
        else:   
            self.backGrounded = False


        #apply the rotation to the wheels
        self.frontRotation += -self.frontAngularVelocity  * dt
        self.backRotation += -self.backAngularVelocity * dt 

    def slipRatio(self):
        #calculate the sped of the car in the given direction
        speed = (self.velocity.magnitude() /2) * self.direction

        #set temp values for slip
        frontSlipRatio = 0
        backSlipRatio = 0
        if speed != 0:
            #calculate slip ratio
            frontSlipRatio = (self.frontAngularVelocity * self.wheelRadius - speed) / abs(speed)
            backSlipRatio = (self.backAngularVelocity * self.wheelRadius -speed) / abs(speed)

        #calculate b value
        bFront = self.stiffness * frontSlipRatio
        bBack = self.stiffness * backSlipRatio

        #calculate the friction coefficient
        frontFrictionCoefficient = self.peak * math.sin(self.shape * math.atan(bFront - self.curvature * (bFront - math.atan(bFront))))
        backFrictionCoefficient = self.peak * math.sin(self.shape * math.atan(bBack - self.curvature * (bBack- math.atan(bBack))))

        #return values
        return [frontFrictionCoefficient, backFrictionCoefficient]

    def weightTransfer(self):
        #calculate a percentage of current length / max length
        frontSuspentionRatio = self.frontSpringLength / self.springTravel
        backSuspentionRatio = self.backSpringLength / self.springTravel

        #find the total travel used between both springs
        total = frontSuspentionRatio + backSuspentionRatio

        #calculate weight on wheel depending on whats grounded
        if self.frontGrounded and not self.backGrounded:
            self.weightFront = (self.mass * 9.81) * (1 -(frontSuspentionRatio / total)) * 2
        elif not self.frontGrounded and self.backGrounded:
            self.weightBack = (self.mass * 9.81) * (1 - (backSuspentionRatio / total)) * 2
        else:
            self.weightFront = (self.mass * 9.81) * (1 -(frontSuspentionRatio / total))
            self.weightBack = (self.mass * 9.81) * (1 - (backSuspentionRatio / total))

    def suspensionForce(self, startPoint, springLength, staticList):
        #the contact point is found
        raycast = self.raycast(startPoint, -self.up, self.springTravel + self.wheelRadius, staticList)
        contactPoint = raycast[0]
        hit = raycast[1]
        #previous length calculated
        lastLength = springLength
        #current length calculated
        springLength = math.sqrt(((startPoint.x - contactPoint.x) **2) + ((startPoint.y - contactPoint.y) **2)) - self.wheelRadius
        #change in length over time calculated
        springVelocity = (lastLength - springLength) / dt


        #The spring force is calculated by finding the difference in length between the springs outstretched length and springs actual length.
        springForce = self.springStiffness * (self.springTravel - springLength)
        #damper force calculated using the springs velocity and stiffness
        damperForce  = self.damperStiffness * springVelocity

        #suspension force found
        suspensionForce = (springForce + damperForce) * self.up

        #return in a list
        return [suspensionForce, springLength, contactPoint, hit]
    
    def raycast(self, startPoint, direction, length, staticList):
        endPoint = startPoint + (direction * length)
        #set default value for max
        maxLength = math.inf

        #set current contact point
        contactPoint = Vector2(0,0)

        for j in range(0, len(staticList)):
            maxSideA = max(self.width + length, self.height + length)
            maxSideB = max(staticList[j].width, staticList[j].height)
            #check if the distance is less than longest sides added
            if self.distanceToObject(staticList[j]) < (max((maxSideA, maxSideB)) / 2):
                verticesB = staticList[j].vertices

                #check each edge in shape B against the edge used in shape A
                for i in range(0, len(verticesB)):
                    #set the first and second point to make an edge in shape B
                    vertexB1 = verticesB[i]
                    vertexB2 = verticesB[(i +1) % len(verticesB)]

                    #compare lines and return true if they intersect 
                    if self.intersectSegments(startPoint, endPoint, vertexB1, vertexB2):
                        #calculate gradients
                        
                        #make sure the gradients are not the same
                        if endPoint.x - startPoint.x == 0 : #prevent divide by 0 error
                            #set x
                            x = startPoint.x
                            #set gradients
                            mB = (vertexB2.y - vertexB1.y) / (vertexB2.x - vertexB1.x)
                            cB = vertexB1.y - (mB * vertexB1.x)
                            #set y
                            y = (mB * x) + cB
                        elif vertexB2.x - vertexB1.x == 0:
                            #set x
                            x = vertexB1.x 
                            #set gradients
                            mA = (endPoint.y - startPoint.y) / (endPoint.x - startPoint.x)
                            cA = startPoint.y - (mA * startPoint.x)
                            #set y
                            y = (mA * x) + cA
                        else:
                            #calculate gradients
                            mA = (endPoint.y - startPoint.y) / (endPoint.x - startPoint.x)
                            mB = (vertexB2.y - vertexB1.y) / (vertexB2.x - vertexB1.x)
                        
                            #calculate y intercepts
                            cA = startPoint.y - (mA * startPoint.x)
                            cB = vertexB1.y - (mB * vertexB1.x)

                            #find the x and y
                            x = (cB - cA) / (mA - mB)
                            y = (mA * x) + cA

                        #calculate length 
                        ab = startPoint - Vector2(x,y)
                        length = math.sqrt((ab.x **2) + (ab.y **2))
                        
                        #find the point closest to the car
                        if length < maxLength:
                            maxLength = length
                            contactPoint = Vector2(x,y)

                        self.groundNormal = (vertexB2 - vertexB1).normalize()

        if contactPoint != Vector2(0,0):
            output = [contactPoint, True]
            return output

        output = [endPoint, False]
        return output
    
    def fuel(self):
        #check if fuel has ran out
        if self.fuelTimer > self.fuelTime:
            #out of fuel
            self.outOfFuel = True
        else:
            #increase the fuel timer
            if self.startGame:
                self.fuelTimer += timePassed / 1000  #increase timer

        #draw fuel gauge image
        screen.blit(self.fuelGuage, (screenWidth/2-107, 740))

        #calculate the angle of the needle using the fuel remaining
        fuelPercentage = self.fuelTimer / self.fuelTime
        needleAngle = fuelPercentage * 180

        needleImg = pygame.transform.rotate(self.fuelNeedle, needleAngle - 90)
        #rotate the needle offset angle
        rotatedOffset = self.fuelNeedleOffset.rotate(-needleAngle + 90)

        #set image to the centre of itself
        needleRect = needleImg.get_rect(center = Vector2(screenWidth/2, 860) - rotatedOffset)
        #draw the needle
        screen.blit(needleImg, needleRect)

    def rpmMeter(self):
        #draw the tachometer
        screen.blit(self.tachometer, (screenWidth/2 -350, 700))

        #calculate the angle of the needle based on the rpm of the engine rpm
        rpmPercentage = self.engineRpm / self.maxRpm
        needleAngle = rpmPercentage * 180

        #rotate the needle
        needleImg = pygame.transform.rotate(self.rpmNeedle, -needleAngle + 90)

        #rotate the needle offset angle
        rotatedOffset = self.rpmNeedleOffset.rotate(needleAngle - 90)

        #set image to the centre of itself
        needleRect = needleImg.get_rect(center = Vector2(screenWidth/2 -257, 800) - rotatedOffset)
        #draw the needle
        screen.blit(needleImg, needleRect)

        #distance text
        distanceRemaining = (70400 - self.position.x) / (self.width / 4)
        
        #rcalculate 
        text = self.rpmFont.render(str(round(distanceRemaining)) + "m", True, (255,201, 14))
        #render the distance remaining
        screen.blit(text, (screenWidth/2-285,810))

    def countdown(self):
        self.countdownTimer -= timePassed / 1000

        #countdown text
        if round(self.countdownTimer) <= 0:
            self.startGame = True
            self.startCountdown = False
        
        text = self.countdownFont.render(str(round(self.countdownTimer)), True, (0,0,0))

        screen.blit(text, (screenWidth/2 -50, screenHeight/2 - 400))

    def timeText(self):
        #incriment timer with time
        self.timer += timePassed / 1000

        #calculate the miliseconds, seconds and minutes
        miliseconds = (self.timer * 100) % 100
        seconds = self.timer % 60
        minutes = int(self.timer /60) % 60
        
        #concaternate the text 
        text = self.timeFont.render(f"{minutes:02}:{round(seconds):02}:{round(miliseconds):02}", True, (255,255,255))

        #render the text 
        screen.blit(text, (screenWidth/2 -70, screenHeight/2 + 240))

    def intersectSegments(self, a,b,c,d):
        #either acd should be counter-clockwise or bcd but not both. same for abc and abd 
        #returns true if two segments intersect
        return self.ccwCheck(a,c,d) != self.ccwCheck(b,c,d) and self.ccwCheck(a,b,c) != self.ccwCheck(a,b,d)

    def ccwCheck(self, a,b,c):
        #check if a, b and c are in a counter clockwise direction. returns true if they are
        return (c.y-a.y)*(b.x-a.x) > (b.y-a.y)*(c.x-a.x)
    
    def inputs(self): #move to apply forces (used for dubugging)
        keys = pygame.key.get_pressed()
        if keys[pygame.K_a]:
            #self.addForce(Vector2(-1000,0))
            self.inputHorizontal = -1
        elif keys[pygame.K_d]:
            #self.addForce(Vector2(1000,0))
            self.inputHorizontal = 1 
        else:
            self.inputHorizontal = 0

        if keys[pygame.K_w]:
            #self.addForce(Vector2(-1000,0))
            self.addForce(Vector2(0,-1000))

        if keys[pygame.K_r]: #for debugging
            self.addTorque(-90000)
        
    def distanceToObject(self, objectB):
        #calculate the vector from a to b
        ab = objectB.position - self.position

        #return the distance 
        return math.sqrt(ab.x **2 + ab.y **2)

    def lerp(self, current, goal, time):
        #calculate the gap 
        difference = goal - current

        #move towards the target 
        if difference > time:
            return current + time
        if difference < -time:
            return current - time
        
        #return the goal if the same
        return goal

class Player():
    def __init__(self, width, height, name):  
        #vectors
        self.position = Vector2(0,0)
        self.currentPosition = Vector2(0,0)
        self.velocity = Vector2(0,0)
        self.rotation = 0
        self.up = Vector2(0,0)

        #wheels
        self.frontWheelPosition = Vector2(0,0)
        self.backWheelPosition = Vector2(0,0)
        self.currentFrontWheelPosition = Vector2(0,0)
        self.currentBackWheelPosition = Vector2(0,0)
        
        #Wheel Variables
        self.wheelRadius = 20
        self.wheelRotation = 0
        
        self.width = width
        self.height = height

        #images
        self.carImg = pygame.image.load(path + "car2.png").convert_alpha()
        self.carImg = pygame.transform.scale(self.carImg, (self.width, self.height+20)) 

        self.frontWheelImg = pygame.image.load(path + "tire2.png").convert_alpha()
        self.frontWheelImg = pygame.transform.scale(self.frontWheelImg, (self.wheelRadius * 2, self.wheelRadius *2))

        self.backWheelImg = pygame.image.load(path + "tire2.png").convert_alpha()
        self.backWheelImg = pygame.transform.scale(self.backWheelImg, (self.wheelRadius * 2, self.wheelRadius *2))

        #lerp
        self.lerpSpeed = 5

        self.name = name
        self.nameFont = pygame.font.SysFont("Bahnschrift", 30, bold = True)
        self.nameText = self.nameFont.render(self.name, True, (255,255,255))
        
    def renderPlayer(self, offset):
        #calculate the distance between current and taregt position
        gap = (self.currentPosition - self.position).magnitude()
        frontWheelGap = (self.currentFrontWheelPosition - self.frontWheelPosition).magnitude()
        backWheelGap = (self.currentBackWheelPosition - self.backWheelPosition).magnitude()

        #lerp the cars current position to the recieved position
        self.currentPosition.x = self.lerp(self.currentPosition.x, self.position.x, dt * self.lerpSpeed * gap **2)
        self.currentPosition.y = self.lerp(self.currentPosition.y, self.position.y, dt * self.lerpSpeed * gap **2)

        #lerp the front and back wheel positions
        self.currentFrontWheelPosition.x = self.lerp(self.currentFrontWheelPosition.x, self.frontWheelPosition.x, dt * self.lerpSpeed * frontWheelGap **2)
        self.currentFrontWheelPosition.y = self.lerp(self.currentFrontWheelPosition.y, self.frontWheelPosition.y, dt * self.lerpSpeed * frontWheelGap **2)
        self.currentBackWheelPosition.x = self.lerp(self.currentBackWheelPosition.x, self.backWheelPosition.x, dt * self.lerpSpeed * backWheelGap **2)
        self.currentBackWheelPosition.y = self.lerp(self.currentBackWheelPosition.y, self.backWheelPosition.y, dt * self.lerpSpeed * backWheelGap **2)

        #rotate the car image with the given rotation
        carImg = pygame.transform.rotate(self.carImg, math.degrees(self.rotation))

        #calculate the wheels speed
        angularVelocity = (self.velocity.x / self.wheelRadius /2)
        self.wheelRotation += angularVelocity * dt

        #rotate the wheels based on the cars speed
        frontWheelImg = pygame.transform.rotate(self.frontWheelImg, math.degrees(self.wheelRotation))
        backWheelImg = pygame.transform.rotate(self.backWheelImg, math.degrees(self.wheelRotation))

        #set the cars rectangle from the image
        carRect = carImg.get_rect(center = self.currentPosition +10 * self.up - offset)
        #render the car image
        screen.blit(carImg, carRect)

        #find centre point of wheel
        frontWheelRect = frontWheelImg.get_rect(center = self.currentFrontWheelPosition - offset  - self.wheelRadius * -self.up)
        backWheelRect = backWheelImg.get_rect(center = self.currentBackWheelPosition - offset - self.wheelRadius * -self.up)

        #draw the wheels
        screen.blit(frontWheelImg, frontWheelRect)
        screen.blit(backWheelImg, backWheelRect)

        screen.blit(self.nameText, self.currentPosition - Vector2(50,60) - offset)

    def lerp(self, current, goal, time):
        #calculate the gap 
        difference = goal - current

        #move towards the target 
        if difference > time:
            return current + time
        if difference < -time:
            return current - time
        
        #return the goal if the same
        return goal

class fuelCan():
    def __init__(self, x, y):
        self.position = Vector2(x,y)
        self.fuelCan = pygame.image.load(path + "fuel-canister.png").convert_alpha()
        self.fuelCanRect = self.fuelCan.get_rect(center = self.position)

    def drawCan(self, offset):
        self.fuelCanRect = self.fuelCan.get_rect(center = self.position - offset)
        screen.blit(self.fuelCan, self.fuelCanRect)

    def collision(self, carRect):
        if self.fuelCanRect.colliderect(carRect):
            return True
        
class Button():
    def __init__(self, x, y, width, height, text, fontSize, bold, colour, highlightedColour, textColour):
        #position variables
        self.x = x
        self.y = y
        
        #diemntions of button
        self.width = width
        self.height = height

        #text colour variables / size
        self.text = text
        self.fontSize = fontSize
        self.colour = (0,0,0)
        self.normalColour = colour
        self.highlightedColour = highlightedColour
        self.textColour = textColour
        
        #font variables
        self.font = pygame.font.SysFont("Bahnschrift", self.fontSize, bold = bold)
        self.fontWidth , self.fontHeight = self.font.size(text)
        self.buttonRect = pygame.Rect(self.x, self.y, self.width, self.height)

    def drawButtonRaw(self):
        #define the button rectange
        self.buttonRect = pygame.Rect(self.x, self.y, self.width, self.height)
        #draw the rectangle
        pygame.draw.rect(screen, self.colour, self.buttonRect)

        #render the text 
        text = self.font.render(self.text, True, self.textColour)
        #center the text in the middle of the button
        screen.blit(text,(self.x +((self.width/2) - (self.fontWidth/2)), self.y + ((self.height/2) - self.fontHeight/2)))
    
    def clickButton(self):
        #get mouse position
        mousePos = pygame.mouse.get_pos()
        #compare mouse position to button rectangle
        if self.buttonRect.collidepoint(mousePos):
            #has collided
            self.colour = self.highlightedColour
            return True
        else:
            #has not collided
            self.colour = self.normalColour

class Text():
    def __init__(self, text, textSize, colour, bold, x, y):
        #set the font and text 
        self.font = pygame.font.SysFont("Bahnschrift", textSize, bold = bold)
        self.text = self.font.render(text, True, colour)

        #get the font width and height
        self.fontWidth, self.fontHeight = self.font.size(text)

        #set the position
        self.x = x
        self.y = y

    def renderText(self):
        #render text at the center of itself
        screen.blit(self.text, (self.x - (self.fontWidth/2), self.y))

class InputBox():
    def __init__(self, x, y, width, height, colour, highlightedColour, fontSize, bold, textColour):
        self.x = x
        self.y = y

        #dimentions
        self.width = width
        self.height = height
        #font
        self.fontSize = fontSize
        self.colour = (0,0,0)
        self.text = ""
        self.normalColour = colour
        self.highlightedColour = highlightedColour
        self.textColour = textColour
        self.font = pygame.font.SysFont("Bahnschrift", self.fontSize, bold = bold)
        self.rect = pygame.Rect(self.x, self.y, self.width, self.height)

        #toggles
        self.toggled = False
        self.untoggleBox()

    def drawBox(self):
        #set the rectangle
        self.rect = pygame.Rect(self.x, self.y, self.width, self.height)
        #draw the rect
        pygame.draw.rect(screen, self.colour, self.rect)

        #get width and height of text 
        fontWidth , fontHeight = self.font.size(self.text)

        #set text
        text = self.font.render(self.text, True, self.textColour)
        #centralise and draw text 
        screen.blit(text,(self.x +((self.width/2) - (fontWidth/2)), self.y + ((self.height/2) - fontHeight/2)))

    def hoverBox(self):
        #get mosue position
        mousePos = pygame.mouse.get_pos()

        #if mouse inside box return true
        if self.rect.collidepoint(mousePos):
            return True
        else:
            return False
    
    def toggleBox(self):
        #set toggle to be true
        self.toggled = True
        #change colour
        self.colour = self.highlightedColour
    
    def untoggleBox(self):
        #toggle false
        self.toggled = False
        #change colour
        self.colour = self.normalColour

class Network():
    def __init__(self):
        #image path
        f = open(pathToSettings, "rt")
        ip = ""
        for i in f.readlines()[2]:
            if i == ";":
                break
            ip += i 
        f.close()

        self.ip = ip #str(input("Enter Sever ip: "))
        self.clientIp = socket.gethostbyname(socket.gethostname())
        self.clientPort = random.randint(7000, 9000)
        self.serverPort = 9999

        self.connected = False

        self.server = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

        self.requests = queue.Queue()

        self.lobbyList = []
        self.lobbyName = ""

        self.playerLookup = {} #dictionary of ip to names   -  ip : name

        self.players = {} #dictionary of ip to player objects ip : player object

        self.startGame = False

    def recieveData(self):
        while True:
            try:
                #get the data and address from any messages recieved
                data, _ = self.server.recvfrom(1024)

                #push the data and address into the queue
                self.requests.put(data.decode())

            except:
                pass
            
            if threadingEvent.is_set():
                break

            time.sleep(0.001)

    def dataCheck(self):
        while True:
            while not self.requests.empty():
                request = self.requests.get()
                requestList = request.split(' ')

                if self.startGame:
                    if request.startswith("RP"):
                        #recieve position and rotation
                        
                        #set the other players ip
                        playerIp = requestList[1]
                        #update player lookup
                        player = self.players[playerIp]

                        #set all sent information
                        player.position = Vector2(float(requestList[2]), float(requestList[3]))
                        player.rotation = float(requestList[4])
                        player.up = Vector2(float(requestList[5]), float(requestList[6]))
                        player.frontWheelPosition = Vector2(float(requestList[7]), float(requestList[8]))
                        player.backWheelPosition = Vector2(float(requestList[9]), float(requestList[10]))
                        player.velocity.x = float(requestList[11])
                    
                    elif request.startswith("STARTCOUNTDOWN"):
                        car.startCountdown = True

                else:

                    if request.startswith("SL"):
                        #recieve lobbys
                        print("recieving lobbys")
                        self.lobbyList = []
                        for i in range(1, len(requestList)-1):
                            self.lobbyList.insert(i, requestList[i])

                    elif request.startswith("SN"):
                        #send names
                        #reset player lookup
                        self.playerLookup = {}

                        #loop though names and update the player lookup
                        for i in range(1, len(requestList)-1,2):
                            self.playerLookup.update({requestList[i] : requestList[i+1]})
                            #put packet totgrather

                    elif request.startswith("JS"):
                        self.connected = True
                        print("Joined server")

                    elif request == "START":
                        self.startGame = True
                        #loop through all player ip addresses 
                        for ip in self.playerLookup.keys():
                            #Make sure not our own ip address
                            if ip != self.clientIp:
                                #create a new player and update the player object dictionary 
                                newPlayer = Player(160, 50, self.playerLookup[ip])
                                self.players.update({ip : newPlayer})
                                print(self.players)
                
                if threadingEvent.is_set():
                    break
            time.sleep(0.001)
        
    def connectToServer(self, name):
        #try to connect to server
        try:
            #send name to server
            self.sendName(name)
        except socket.error as e:
            #cant join server
            print("Unable to join server. Error: ", e)
    
    def createLobby(self, lobbyName, lobbySize):
        #make sure there is a lobby name
        if lobbyName == "":
            self.lobbyName = "Lobby" + str(random.randint(1000,9999))
        else:
            self.lobbyName = lobbyName

        #send lobby name and size
        data = "CL " + str(lobbySize) + " " + str(self.lobbyName)
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def joinLobby(self):
        #concaternate data
        data = "JL " + str(self.lobbyName)
        #send to server
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def getLobbyList(self):
        data = "RL" 
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def sendName(self, name):
        #make sure name is not empty 
        if name == "":
            name = "Player" + str(random.randint(1000,9000))

        #concaternate data
        data = "JS " + name 
        #send name to server
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def leaveLobby(self):
        data = "LL " + self.lobbyName #+ " " + self.clientIp + " " + str(self.clientPort)
        print("leaving lobby")
        self.server.sendto(data.encode(), (self.ip, self.serverPort))
        self.lobbyName = ""

    def readyUp(self):
        data = "READY " 
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def gameReady(self):
        data = "LG "
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

    def sendPlayerInformation(self, position, rotation, up, frontWheel, backWheel, xVelocity):
        #concaternate information
        data = "RP " + str(position.x) + " " + str(position.y) + " " + str(rotation) +  " " + str(up.x) + " " + str(up.y) +  " " + str(frontWheel.x) + " " + str(frontWheel.y)+  " " + str(backWheel.x) + " " + str(backWheel.y) + " " + str(xVelocity)
        #send information
        self.server.sendto(data.encode(), (self.ip, self.serverPort))

#create the network object 
global networkManager
networkManager = Network()


#-----------------Menus------------------
def mainMenu():
    #create button objects
    playButton = Button(screenWidth/2 - 50, screenHeight/2 - 70, 100, 70, "Play", 30, False, (0,255,0),(0,200,0), (0,0,0))
    connectButton = Button(screenWidth/2 - 50, screenHeight/2 - 70, 100, 70, "Connect", 30, False, (0,255,0),(0,200,0), (0,0,0))
    exitButton = Button(screenWidth/2 - 50, screenHeight/2 - 70 + 80, 100, 70, "Exit", 30, False, (255,0,0), (200,0,0),(0,0,0))
    
    #Create text objects
    titleText = Text("Racing Game", 60, (0,0,0), True, screenWidth/2, 100)
    
    
    enterNameText = Text("Enter your name", 30, (0,0,0), False, screenWidth/2, 230)

    nameInputBox = InputBox(550, 275, 500, 50, (0,38,77), (0,26,51), 30, False, (255,255,255))

    menu = True
    while menu:
        #Make it so the game window can close
        click = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                click =  True
                nameInputBox.untoggleBox()
            if event.type == pygame.TEXTINPUT and nameInputBox.toggled:
                #make sure you cannot enter space
                if event.text != " ":
                    nameInputBox.text += event.text
            if event.type == pygame.KEYDOWN and nameInputBox.toggled:
                if event.key == pygame.K_BACKSPACE:
                    nameInputBox.text = nameInputBox.text[:-1]


        if not menu:
            break
        
        #fill screen
        screen.fill((200, 200, 200))

        titleText.renderText()
        
        enterNameText.renderText()


        if not networkManager.connected:
            connectButton.drawButtonRaw()
            if connectButton.clickButton() and click:
                networkManager.connectToServer(nameInputBox.text)
        else:
            playButton.drawButtonRaw()
            if playButton.clickButton()and click:
                playMenu()
                menu = False
                break
        
        if nameInputBox.hoverBox() and click:
            nameInputBox.toggleBox()

        if exitButton.clickButton()and click:
            menu = False
            GameManager.exitGame()
        
        nameInputBox.drawBox()
        exitButton.drawButtonRaw()
        


        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def playMenu():
    #create button objects
    garageButton = Button(screenWidth/2 - 50, screenHeight/2 - 70 - 80, 100, 70, "Garage", 30, False, (0,255,0),(0,200,0), (0,0,0))
    createLobbyButton = Button(screenWidth/2 - 90, screenHeight/2 - 70, 180, 70, "Create Lobby", 30, False, (0,255,0),(0,200,0), (0,0,0))
    joinLobbyButton = Button(screenWidth/2 - 75, screenHeight/2 - 70 + 80, 150, 70, "Join Lobby", 30, False, (255,0,0), (200,0,0),(0,0,0))
    
    backButton = Button(100, 800, 100, 70, "Back", 30, False, (255,0,0), (200,0,0),(0,0,0))

    #game loop
    menu = True
    while menu:
        #Make it so the game window can close
        click = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                #exit game
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                #has clicked
                click =  True
        
        #close menu if false
        if not menu:
            break
            
        #fill screen
        screen.fill((200,200,200))

        #check button
        if garageButton.clickButton() and click:
            garageMenu()
            break
        #check button
        if createLobbyButton.clickButton()and click:
            createLobbyMenu()
            break
        #check button
        if joinLobbyButton.clickButton()and click:
            joinLobbyMenu()
            break
        #check button
        if backButton.clickButton() and click:
            mainMenu()
            break
        
        #draw buttons
        backButton.drawButtonRaw()
        garageButton.drawButtonRaw()
        createLobbyButton.drawButtonRaw()
        joinLobbyButton.drawButtonRaw()

        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def joinLobbyMenu(): 
    #get the lobby lsit automatically
    networkManager.getLobbyList()
    menu = True
    #set the text 
    titleText = Text("Join a Lobby", 60, (0,0,0), False, screenWidth/2, 100)
    
    #set the buttons
    refreshButton = Button(800, 800, 100, 70, "Refresh", 30, False, (255,0,0), (200,0,0),(0,0,0))
    backButton = Button(100, 800, 100, 70, "Back", 30, False, (255,0,0), (200,0,0),(0,0,0))

    #load the lobbies
    def loadLobbys():
        buttonLobbyList = []
        lobbyList = networkManager.lobbyList

        #loop though each lobby and create a button for it
        for i in range(0,len(lobbyList)):
            button = Button(screenWidth/2-250/2, 200 + (70 * i), 250, 50, lobbyList[i], 30, False, (255,0,0), (200,0,0),(0,0,0))
            #append button to button list
            buttonLobbyList.append(button)
        #output button list
        return buttonLobbyList
    
    while menu:
        #Make it so the game window can close
        click = False
        #check for events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                click =  True
        
        if not menu:
            break
        #fill screen
        screen.fill((200,200,200))

        #loop though each button and display it
        for lobbyButton in loadLobbys():
            #check if button clicked
            if lobbyButton.clickButton() and click:
                #join a lobby if clicked
                networkManager.lobbyName = lobbyButton.text
                networkManager.joinLobby()
                lobbyMenu()
                break 
            if menu:
                #draw button
                lobbyButton.drawButtonRaw()

        #render text 
        titleText.renderText()

        #check for button click
        if backButton.clickButton() and click:
            playMenu()
            break
        #check for button click
        if refreshButton.clickButton() and click:
            networkManager.getLobbyList()

        #render buttons
        backButton.drawButtonRaw()
        refreshButton.drawButtonRaw()

        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def garageMenu():
    menu = True
    titleText = Text("Garage", 60, (0,0,0), False, screenWidth/2, 100)
    
    backButton = Button(100, 800, 100, 70, "Back", 30, False, (255,0,0), (200,0,0),(0,0,0))
    while menu:
        #Make it so the game window can close
        click = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                click =  True

        if not menu:
            break
        
        screen.fill((100,100,100))

        titleText.renderText()

        if backButton.clickButton() and click:
            playMenu()
            break 

        backButton.drawButtonRaw()

        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def createLobbyMenu():
    menu = True
    #create text
    titleText = Text("Create Lobby", 60, (0,0,0), False, screenWidth/2, 100)
    enterNameText = Text("Enter Lobby Name", 30, (0,0,0), False, screenWidth/2, 230)

    #create input box
    lobbyNameInputBox = InputBox(550, 275, 500, 50, (0,38,77), (0,26,51), 30, False, (255,255,255))
    
    #create buttons
    createLobbyButton = Button(screenWidth/2, screenHeight/2, 200, 70, "Create Lobby", 30, False, (255,0,0), (200,0,0),(0,0,0))
    backButton = Button(100, 800, 100, 70, "Back", 30, False, (255,0,0), (200,0,0),(0,0,0))

    while menu:
        #Make it so the game window can close
        click = False
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                click =  True
                lobbyNameInputBox.untoggleBox()
            if event.type == pygame.TEXTINPUT and lobbyNameInputBox.toggled:
                lobbyNameInputBox.text += event.text
            if event.type == pygame.KEYDOWN and lobbyNameInputBox.toggled:
                if event.key == pygame.K_BACKSPACE:
                    lobbyNameInputBox.text = lobbyNameInputBox.text[:-1]

        if not menu:
            break
    
        #fill  screen
        screen.fill((200,200,200))

        #redner text 
        titleText.renderText()
        enterNameText.renderText()

        #check buttons for click
        if backButton.clickButton() and click:
            playMenu()
            break
        #check buttons for click
        if createLobbyButton.clickButton() and click:
            networkManager.createLobby(lobbyNameInputBox.text, 5)
            lobbyMenu()
            break 
        #check input box for click
        if lobbyNameInputBox.hoverBox() and click:
            lobbyNameInputBox.toggleBox()

        #render buttons and input box
        createLobbyButton.drawButtonRaw()
        lobbyNameInputBox.drawBox()
        backButton.drawButtonRaw()

        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def lobbyMenu():
    menu = True
    #create text 
    titleText = Text(networkManager.lobbyName, 60, (0,0,0), False, screenWidth/2, 100)
    
    #create buttons
    leaveLobbyButton = Button(100, 800, 200, 70, "Leave Lobby", 30, False, (255,0,0), (200,0,0),(0,0,0))
    readyButton = Button(600, 800, 200, 70, "Ready", 30, False, (255,0,0), (200,0,0),(0,0,0))

    def loadNames():
        #create names 
        textNameList = []
        nameList = networkManager.playerLookup
        i =0
        for name in nameList.values():
            #loop through each name in dictionary and add it to the name list 
            nameText = Text(str(name), 30, (0,0,0), False, screenWidth/2, 400 + (30 * i))
            textNameList.append(nameText)
            i+=1
        
        #output name list
        return textNameList
    
    while menu:
        #Make it so the game window can close
        click = False
        #check for events
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                menu = False
                GameManager.exitGame()
            if event.type == pygame.MOUSEBUTTONUP:
                click =  True

        if not menu:
            break
        
        #fill screen
        screen.fill((200,200,200))

        #render title text 
        titleText.renderText()

        #draw all player names
        for nameText in loadNames():
            nameText.renderText()

        #check for button clicked
        if leaveLobbyButton.clickButton() and click:
            networkManager.leaveLobby()
            playMenu()
            break 
        
        #check for button clicked
        if readyButton.clickButton() and click:
            networkManager.readyUp()

        #start game
        if networkManager.startGame:
            print("starting game")
            gameMenu()
            break
        
        #draw buttons
        leaveLobbyButton.drawButtonRaw()
        readyButton.drawButtonRaw()

        #Update display and tick the clock
        pygame.display.update()
        clock.tick(fps)

def gameMenu():
    #Other Settings
    counter = 0
    frameStart = time.time()  

    #create objects    
    global car       
    car = Car(1425,700,0,160,50,1500,(255,0,0))
    gameManager = GameManager(car)

    def start(): #called only once at the start of the game
        gameManager.spawnTerrain()
        networkManager.gameReady()

    def update(): #This subroutine is called every frame
        gameManager.renderObjects()
        
    def fixedUpdate(): #This subroutine is called a fixed number of times a second
        gameManager.updateObjects()

    #call the start method
    start()

    #loop condition
    running = True

    global timePassed
    timePassed = 0

    #game loop
    while running: 
        #set the current time each frame
        currentTime = time.time()

        #increase counter and update the frame start
        counter += currentTime - frameStart
        frameStart = currentTime

        #Counter limit
        if counter > 0.2:
            counter = 0.2

        #Fixed GameLoop
        while counter > dt:
            fixedUpdate()
            counter -= dt

        #normal
        update()

        #Make it so the game window can close
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
                GameManager.exitGame()
        
        if not running:
            break

        #Update display and tick the clock
        pygame.display.update()
        timePassed = clock.tick(fps)

def testMenu():
    #Other Settings
    counter = 0
    frameStart = time.time()  

    testObject = physicsObject(100,500,0,100,100,50,(255,255,255))

    def update(): #This subroutine is called every frame
        screen.fill((50,50,70))
        testObject.renderObject()
        
    def fixedUpdate(): #This subroutine is called a fixed number of times a second
        testObject.updatePhysics()

        testObject.addForce(Vector2(10,0))

    #loop condition
    running = True

    global timePassed
    timePassed = 0

    #game loop
    while running: 
        #set the current time each frame
        currentTime = time.time()

        #increase counter and update the frame start
        counter += currentTime - frameStart
        frameStart = currentTime

        #Counter limit
        if counter > 0.2:
            counter = 0.2

        #Fixed GameLoop
        while counter > dt:
            fixedUpdate()
            counter -= dt

        #normal
        update()

        #Make it so the game window can close
        for event in pygame.event.get():
            if event.type == pygame.QUIT:
                running = False
                GameManager.exitGame()
        
        if not running:
            break

        #Update display and tick the clock
        pygame.display.update()
        timePassed = clock.tick(fps)
threadingEvent = threading.Event()
#start the main menu

#create the threading and start the threads
t1 = threading.Thread(target = networkManager.recieveData)
#start daemon
t1.daemon = True
t1.start()

t2 = threading.Thread(target = networkManager.dataCheck)
#start daemon
t2.daemon = True
t2.start()

mainMenu()
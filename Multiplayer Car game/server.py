import socket
import threading
import queue
import time

ip = socket.gethostbyname(socket.gethostname())
port = 9999

server = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)

#set requests queue
requests = queue.Queue()

lobbys = []
lobbyLookup = {}
names = {}

class Lobby:
    def __init__(self, lobbyName, lobbySize):
        #player list 
        self.players = []

        #ready and start dictionary
        self.readyDict = {}
        self.startDict = {}
        #name of lobby 
        self.lobbyName = lobbyName
        self.maxLobbySize = lobbySize
        self.currentLobbySize = 1

        #check if game has started
        self.gameStarted = False

    def addPlayer(self, player):
        #add player to the list 
        self.players.append(player)
        #incriment lobby size
        self.currentLobbySize += 1

    def removePlayer(self, player):
        #loop though each player until found
        for i in self.players:
            if i == player:
                #remove player
                self.players.remove(i)

#try to bind the server with the given ip and port
try:
    server.bind((ip, port))
    print("Server Started at " + ip)
except socket.error as e:
    print(str(e))

#recieves data and adds it to a queue
def recieve():
    while True:
        try:
            #get the data and address from any messages recieved
            data, addr = server.recvfrom(1024)
            request = data.decode()
            #push the data and address into the queue
            requests.put((request, addr))
        except:
            pass
        #slow down ever so slightly to reduce cpu usage
        time.sleep(0.001)

def sendPlayerNames(lobby):
    #create a list of player names
    data = "SN "
    for player in lobby.players:
        #send names
        data += str(player[0]) + " " + names[player[0]] + " " 
                        
    #send the list of player names to each player
    for player in lobby.players:
        sendTo(data, player)

def checkReady(lobby):
    for value in lobby.readyDict.values():
        if not value:
            return False
    return True

def checkAllPlayersIn(lobby):
    for value in lobby.startDict.values():
        if not value:
            return False
    return True

def dataCheck():
    while True:
        #check is request queue is not empty 
        while not requests.empty():
            #dequeue the list and get the request
            request, addr = requests.get()

            #split the reqeust up into a list at each space
            requestList = request.split(' ')

            if request.startswith("RP"):
                #recieve position
                lobby = lobbyLookup[addr[0]]

                #send ip address
                data = "RP " + str(addr[0]) + " "

                #concaternate data again
                for i in range(1,len(requestList)):
                    data += requestList[i] + " " 
                
                #send data to all players
                for player in lobby.players:
                    #make sure not to send to yourself
                    if player != addr:
                        sendTo(data, player)

            elif request.startswith("CL"):
                #create lobby
                #get the lobby name
                lobbyName = requestList[2]
                #print to console
                print(names[addr[0]] + " made a new lobby called " + lobbyName)

                #create a new lobby 
                lobby = Lobby(lobbyName, requestList[1])

                #append the lobby 
                lobbys.append(lobby)

                #add player to lobby 
                lobby.addPlayer(addr)
                #update the lobby lookup
                lobbyLookup.update({addr[0] : lobby}) 

                #update the ready dictionary to false
                lobby.readyDict.update({addr[0] : False})
                lobby.startDict.update({addr[0] : False})
                
                #print to console
                print(names[addr[0]] + " has joined " + lobbyName)

                #send all players in the lobby a new name list 
                sendPlayerNames(lobby)


            elif request.startswith("JL"):
                #join lobby
                #find the lobby 
                for lobby in lobbys:
                    if requestList[1] == lobby.lobbyName:
                        #add player to lobby
                        lobby.addPlayer(addr)
                        #update lobby look up
                        lobbyLookup.update({addr[0] : lobby}) 

                        #update the ready dictionary to false
                        lobby.readyDict.update({addr[0] : False})
                        lobby.startDict.update({addr[0] : False})

                        #print to console
                        print(names[addr[0]] + " has joined " + requestList[1])

                        #send new list of players in lobby
                        sendPlayerNames(lobby)
                
            elif request == "RL":
                #request lobbys
                data = "SL "
                for lobby in lobbys:
                    if not lobby.gameStarted:
                        data += lobby.lobbyName + " "
                sendTo(str(data), addr)
            
            elif request.startswith("LL"):
                #leave lobby
                i = 0
                for lobby in lobbys:
                    if requestList[1] == lobby.lobbyName:
                        lobby.removePlayer(addr)
                        del lobbyLookup[addr[0]]
                        print(names[addr[0]] + " has left", requestList[1])
                        #if not lobbyLookup:
                           # print(requestList[1] + "  is empty. Lobby will be removed.")
                            #del lobbys[i]
                        sendPlayerNames(lobby)
                    i+=1

            elif request.startswith("JS"):
                #send back a join server request to client
                sendTo("JS", addr)
                #update names lookup dictionary 
                names.update({addr[0] : requestList[1]}) 

                print(names[addr[0]] + " has connected to the server")

            elif request.startswith("READY"):
                lobby = lobbyLookup[addr[0]]
                lobby.readyDict.update({addr[0] : True})
                print("Starting game in lobby " + lobby.lobbyName)
                print(lobby.players)
                if checkReady(lobby):
                    for player in lobby.players:
                        sendTo("START", player)
                    lobby.gameStarted = True
            
            elif request.startswith("LG"):
                #loaded game
                lobby = lobbyLookup[addr[0]]
                lobby.startDict.update({addr[0] : True})
                print(lobby.players)
                if checkAllPlayersIn(lobby):
                    for player in lobby.players:
                        sendTo("STARTCOUNTDOWN", player)


        time.sleep(0.001)

def sendTo(data, addr):
    try:
        server.sendto(data.encode(), addr)
    except:
        pass

#create thread 1 - recieve data
t1 = threading.Thread(target = recieve)
t1.start()

#create thread 2 - check data
t2 = threading.Thread(target = dataCheck)
t2.start()
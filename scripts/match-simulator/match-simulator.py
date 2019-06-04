## Python Stochastic Simulator Example ##
## Run python scripts using Python 3.5 or newer ##

''' Simulator script:
    ----------------
    Object-Oriented approach to simulating competitive match results.
    
    Author:     Austin Ladshaw
    Date:       06/03/2019
    Copyright:  This software was designed and built by Austin Ladshaw.
                Copyright (c) 2019, all rights reserved.
'''

import random

class MTGData(object):
    def __init__(self):
        self.match_ups = {}
        self.meta_share = {}

    def regDecks(self, list):
        for deck in list:
            self.match_ups[deck] = {}
            self.meta_share[deck] = 0.0

    def readMatchupRates(self, file):
        data = open(file,"r")
        deck_list = []
        deck_set = False
        i = 0
        for line in data:
            strings = line.split("\t")
            #First line is header, so we skip over it
            if i == 0:
                if strings[0].split("\n")[0] != "WinRates":
                    raise ValueError("Unanticipated file format!")
            #Second line list all decks
            elif i == 1:
                if deck_set == False:
                    deck_set = True
                    j = 0
                    for word in strings:
                        if j != 0:
                            deck_list.append(word)
                        j += 1
                    deck_list[len(strings)-2] = deck_list[len(strings)-2].split("\n")[0]
                    self.regDecks(deck_list)
            else:
                deck1 = strings[0]
                if deck1 != "" and deck1 != "\n":
                    for k in range(1,len(strings)-1):
                        self.set_MatchUp(deck1,deck_list[k-1],float(strings[k]))
                    self.set_MatchUp(deck1,deck_list[len(strings)-2],float(strings[len(strings)-1].split("\n")[0]))
            i += 1
        data.close()
    
    def readMetaShare(self, file):
        data = open(file, "r")
        i = 0
        for line in data:
            strings = line.split("\t")
            if i == 0:
                if strings[0].split("\n")[0] != "Meta Game":
                    raise ValueError("Unanticipated file format!")
            else:
                self.set_MetaShare(strings[0],float(strings[1].split("\n")[0]))
            i += 1
        
        data.close()
    
    def set_MatchUp(self, deck1, deck2, d1_v_d2_rate):
        self.match_ups[deck1][deck2] = d1_v_d2_rate
        self.match_ups[deck2][deck1] = 1.0-d1_v_d2_rate

    def set_MetaShare(self, deck, share):
        self.meta_share[deck] = share

    def get_Probability(self, deck1, deck2):
        try:
            return self.match_ups[deck1][deck2]
        except:
            return 0.5

    def get_MetaShare(self, deck):
        try:
            return self.meta_share[deck]
        except:
            return 0.001

    def get_MetaMap(self):
        return self.meta_share

    def __str__(self):
        string = "Meta Game:\n"
        string+= "----------\n"
        for deck in self.meta_share:
            string+= deck + "\t:\t" + str(self.meta_share[deck]) + "\n"
        string += "\nMatch Ups:\n"
        string +=   "----------\n(Decks)"
        first_line = True
        for deck1 in self.match_ups:
            if first_line == True:
                for deck2 in self.match_ups[deck1]:
                    string += "\t" + deck2
                string += "\n"
            string += deck1
            for deck2 in self.match_ups[deck1]:
                string += "\t" + str(self.match_ups[deck1][deck2])
            first_line = False
            string += "\n"
        return string

class Player(object):
    def __init__(self, name, deck):
        # Assumes name and deck are strings
        self.name = name
        self.deck = deck
        self.m_wins = 0
        self.m_losses = 0
        self.m_draws = 0
        
    def get_name(self):
        return self.name
    
    def get_deck(self):
        return self.deck
    
    def get_record(self):
        return (self.m_wins, self.m_losses, self.m_draws)
    
    def update_record(self, win, loss, draw):
        self.m_wins += win
        self.m_losses += loss
        self.m_draws += draw
    
    def erase_record(self):
        self.m_wins = 0
        self.m_losses = 0
        self.m_draws = 0
    
    #Assumes data is Modern Data already filled out
    def simulate_game(self, other, data):
        result = random.uniform(0.0, 1.0)
        prob = data.get_Probability(self.deck, other.deck)
        if result > prob:
            self.m_losses += 1
            other.update_record(1,0,0)
        elif result < prob:
            self.m_wins += 1
            other.update_record(0,1,0)
        else:
            self.m_draws += 1
            other.update_record(0,0,1)

    def __str__(self):
        return (self.name + " : " + self.deck + "\tTourney Record: " + str(self.m_wins) + "-" + str(self.m_losses) + "-" + str(self.m_draws))

class Tourney(object):
    def __init__(self):
        self.players = {}
        self.round = 0
        self.max_rounds = 0
        self.match = {}  #Holds the match data for each round --> [player1, player2, winner]
        #random.seed(0)
    
    def __str__(self):
        string = ""
        string += "\n"
        for round in self.match:
            string += "Round " + str(round) + ":\n"
            for pair in self.match[round]:
                string += "\t" + str(pair[0]) + " v " + str(pair[1]) + " ==> Winner = " + str(pair[2]) + "\n"
            string += "\n"
        for player in self.players:
            string += str(self.players[player]) + "\n"
        return string
    
    def get_player(self, name):
        return self.players[name]
    
    def get_player_map(self):
        return self.players

    def regPlayer(self, player):
        #Assumes a Player object is passed as argument
        self.players[player.get_name()] = player
    
    def regRandomPlayers(self, num, data):
        #Registers a number of random players based on meta data
        sum = 0
        range_map = {}
        for deck in data.get_MetaMap():
            range_map[deck] = []
            range_map[deck].append(sum)
            sum += data.get_MetaShare(deck)
            range_map[deck].append(sum)
        for n in range(0, num):
            name = "R-" + str(n).zfill(len(str(num)))
            rand = random.uniform(0.0, sum)
            res = ""
            for d in range_map:
                if rand >= range_map[d][0] and rand < range_map[d][1]:
                    res = d
                    break
            self.players[name] = Player(name, res)
    
    def determine_rounds(self):
        count = len(self.players)
        if count <= 4:
            self.max_rounds = 2
        elif count <= 8:
            self.max_rounds = 3
        elif count <= 16:
            self.max_rounds = 4
        elif count <= 32:
            self.max_rounds = 5
        elif count <= 64:
            self.max_rounds = 6
        elif count <= 128:
            self.max_rounds = 7
        elif count <= 226:
            self.max_rounds = 8
        elif count <= 409:
            self.max_rounds = 9
        else:
            self.max_rounds = 10

    def RoundPairings(self, round):
        self.round = round
        self.match[round] = []
        list = []
        odd = []

        #Prepare lists for random pairs
        for n in reversed(range(0,round)):
            temp_list = []
            for player in self.players:
                if self.players[player].get_record()[0] == n:
                    temp_list.append(player)
            list.append(temp_list)

        #Choose random pairs
        i = 0
        for sub in list:
            odd.append(False)
            while sub:
                try:
                    player1 = random.choice(sub)
                    sub.remove(player1)
                except:
                    raise IndexError('Player set empty!')
                try:
                    player2 = random.choice(sub)
                    sub.remove(player2)
                    self.match[round].append([player1,player2,None])
                except:
                    odd[i] = True
            #End While loop
            i += 1
            if odd[i-1] == True:
                if len(list) > i:
                    try:
                        player2 = random.choice(list[i])
                        list[i].remove(player2)
                        self.match[round].append([player1,player2,None])
                    except:
                        self.match[round].append([player1,None,None])
                else:
                    self.match[round].append([player1,None,None])

    def simulate_round(self, round, data):
        self.RoundPairings(round)
        for pair in self.match[round]:
            if pair[0] == None or pair[1] == None:
                if pair[0] != None:
                    self.players[pair[0]].update_record(1,0,0)
                    pair[2] = pair[0]
                else:
                    self.players[pair[1]].update_record(1,0,0)
                    pair[2] = pair[1]
            else:
                p1_w = self.players[pair[0]].get_record()[0]
                p2_w = self.players[pair[1]].get_record()[0]
                self.players[pair[0]].simulate_game(self.players[pair[1]], data)
                if p1_w < self.players[pair[0]].get_record()[0]:
                    pair[2] = pair[0]
                if p2_w < self.players[pair[1]].get_record()[0]:
                    pair[2] = pair[1]

    def simulate_tournament(self, data):
        self.determine_rounds()
        for n in range(1, self.max_rounds+1):
            self.simulate_round(n,data)

class TourneyOutcomeSimulator(object):
    def __init__(self, events, data):
        self.tourney = [Tourney()]*events
        self.tot_events = events
        self.player_record = {}
        self.deck_record = {}
        self.mtg_data = data
    
    def __str__(self):
        string = "Results after " + str(self.tot_events) + " tournament(s) simulated of " + str(self.tourney[0].max_rounds) + " rounds each...\n\n"
        i = 0
        for name in self.player_record:
            if i == 0:
                string += "Player\tWins\tLosses\tDraws\tWinRate(%)\n"
            string += name + "(" + self.tourney[0].get_player(name).get_deck() + ")" + "\t" + str(self.player_record[name][0]) + "\t" + str(self.player_record[name][1]) + "\t" + str(self.player_record[name][2]) + "\t" + str(float(self.player_record[name][0])/float(self.player_record[name][1]+self.player_record[name][0]+self.player_record[name][2])*100.0) + "\n"
            i += 1
        string += "\n"
        i = 0
        for deck in self.deck_record:
            if i == 0:
                string += "Deck\tWins\tLosses\tDraws\tWinRate(%)\n"
            string += deck + "\t" + str(self.deck_record[deck][0]) + "\t" + str(self.deck_record[deck][1]) + "\t" + str(self.deck_record[deck][2]) + "\t" + str(float(self.deck_record[deck][0])/float(self.deck_record[deck][1]+self.deck_record[deck][0]+self.deck_record[deck][2])*100.0) + "\n"
            i += 1
        return string

    def regPlayer(self, player):
        self.player_record[player.get_name()] = [0,0,0]
        self.deck_record[player.get_deck()] = [0,0,0]
        for n in range(0,len(self.tourney)):
            self.tourney[n].regPlayer(player)

    def regRandomPlayers(self, num):
        #Registers a number of random players based on meta data
        sum = 0
        range_map = {}
        for deck in self.mtg_data.get_MetaMap():
            range_map[deck] = []
            range_map[deck].append(sum)
            sum += self.mtg_data.get_MetaShare(deck)
            range_map[deck].append(sum)
        for n in range(0, num):
            name = "R-" + str(n).zfill(len(str(num)))
            rand = random.uniform(0.0, sum)
            res = ""
            for d in range_map:
                if rand >= range_map[d][0] and rand < range_map[d][1]:
                    res = d
                    break
        
            self.player_record[name] = [0,0,0]
            self.deck_record[res] = [0,0,0]
            for m in range(0,len(self.tourney)):
                self.tourney[m].regPlayer(Player(name, res))

    def run_simulation(self):
        for n in range(0,len(self.tourney)):
            self.tourney[n].simulate_tournament(self.mtg_data)
            for name in self.player_record:
                self.player_record[name][0] += self.tourney[n].get_player(name).get_record()[0]
                self.player_record[name][1] += self.tourney[n].get_player(name).get_record()[1]
                self.player_record[name][2] += self.tourney[n].get_player(name).get_record()[2]
                self.deck_record[self.tourney[n].get_player(name).get_deck()][0] += self.tourney[n].get_player(name).get_record()[0]
                self.deck_record[self.tourney[n].get_player(name).get_deck()][1] += self.tourney[n].get_player(name).get_record()[1]
                self.deck_record[self.tourney[n].get_player(name).get_deck()][2] += self.tourney[n].get_player(name).get_record()[2]
                self.tourney[n].get_player(name).erase_record()

class RandomTourneyOutcomeSimulator(object):
    def __init__(self, tourneys, events_per_tourney, tot_players, data):
        self.mtg_data = data
        self.rand_tour = [TourneyOutcomeSimulator(events_per_tourney, data)]*tourneys
        self.tot_sims = tourneys
        self.player_record = {}  #Only for static players
        self.deck_record = {}    #For all decks (static and random)
        self.static_players = []
        self.num_players = tot_players

    def __str__(self):
        string = "Results after " + str(self.tot_sims) + " simulations(s) of " + str(self.rand_tour[0].tot_events) + " tournament(s) of " + str(self.rand_tour[0].tourney[0].max_rounds) + " round(s) each...\n\n"
        i = 0
        for name in self.player_record:
            if i == 0:
                string += "Player\tWins\tLosses\tDraws\tWinRate(%)\n"
                string += name + "(" + self.rand_tour[0].tourney[0].get_player(name).get_deck() + ")" + "\t" + str(self.player_record[name][0]) + "\t" + str(self.player_record[name][1]) + "\t" + str(self.player_record[name][2]) + "\t" + str(float(self.player_record[name][0])/float(self.player_record[name][1]+self.player_record[name][0]+self.player_record[name][2])*100.0) + "\n"
            i += 1
        string += "\n"
        i = 0
        for deck in self.deck_record:
            if i == 0:
                string += "Deck\tWins\tLosses\tDraws\tWinRate(%)\n"
                string += deck + "\t" + str(self.deck_record[deck][0]) + "\t" + str(self.deck_record[deck][1]) + "\t" + str(self.deck_record[deck][2]) + "\t" + str(float(self.deck_record[deck][0])/float(self.deck_record[deck][1]+self.deck_record[deck][0]+self.deck_record[deck][2])*100.0) + "\n"
            i += 1
        return string

    def regStaticPlayer(self, player):
        self.static_players.append(player)
        self.player_record[player.get_name()] = [0,0,0]
        self.deck_record[player.get_deck()] = [0,0,0]
        for m in range(0,len(self.rand_tour)):
            for n in range(0,len(self.rand_tour[m].tourney)):
                self.rand_tour[m].tourney[n].regPlayer(player)

    def run_simulations(self):
        num_rand = self.num_players - len(self.static_players)
        for tourney in self.rand_tour:
            for player in self.static_players:
                tourney.regPlayer(player)
            tourney.regRandomPlayers(num_rand)
            
            #Loop through all players for this tourney and put their deck in the deck_record dict if the deck key doesn't yet exist
            
            tourney.run_simulation()
            print(tourney)
            tourney.player_record.clear()
            tourney.deck_record.clear()


## Testing ##
'''
decks = ["G Tron", "UW Control", "Titan", "Dredge", "Jund", "Burn", "Storm"]
data = MTGData()
data.regDecks(decks)
data.set_MetaShare("UW Control", 0.0697)
data.set_MetaShare("G Tron", 0.0659)
data.set_MetaShare("Titan", 0.0386)
data.set_MetaShare("Dredge", 0.0597)
data.set_MetaShare("Jund", 0.0149)
data.set_MetaShare("Burn", 0.0336)
data.set_MetaShare("Storm", 0.0149)

data.set_MatchUp("G Tron", "G Tron", 0.5)
data.set_MatchUp("G Tron", "UW Control", 0.5)
data.set_MatchUp("G Tron", "Titan", 0.405)
data.set_MatchUp("G Tron", "Dredge", 0.529)
data.set_MatchUp("G Tron", "Jund", 0.511)
data.set_MatchUp("G Tron", "Burn", 0.413)
data.set_MatchUp("G Tron", "Storm", 0.368)

data.set_MatchUp("UW Control", "UW Control", 0.5)
data.set_MatchUp("UW Control", "Titan", 0.4)
data.set_MatchUp("UW Control", "Dredge", 0.267)
data.set_MatchUp("UW Control", "Jund", 0.6)
data.set_MatchUp("UW Control", "Burn", 0.559)
data.set_MatchUp("UW Control", "Storm", 0.268)

data.set_MatchUp("Titan", "Titan", 0.5)
data.set_MatchUp("Titan", "Dredge", 0.467)
data.set_MatchUp("Titan", "Jund", 0.525)
data.set_MatchUp("Titan", "Burn", 0.679)
data.set_MatchUp("Titan", "Storm", 0.231)

data.set_MatchUp("Dredge", "Dredge", 0.5)
data.set_MatchUp("Dredge", "Jund", 0.389)
data.set_MatchUp("Dredge", "Burn", 0.652)
data.set_MatchUp("Dredge", "Storm", 0.333)

data.set_MatchUp("Jund", "Jund", 0.5)
data.set_MatchUp("Jund", "Burn", 0.588)
data.set_MatchUp("Jund", "Storm", 0.565)

data.set_MatchUp("Burn", "Burn", 0.5)
data.set_MatchUp("Burn", "Storm", 0.745)

data.set_MatchUp("Storm", "Storm", 0.5)

tourney = Tourney()
tourney.regPlayer(Player("Austin","G Tron"))
tourney.regPlayer(Player("Stephen","UW Control"))
tourney.regPlayer(Player("Rob","Titan"))
tourney.regPlayer(Player("Jon","Dredge"))
tourney.regPlayer(Player("Ken","Jund"))
tourney.regPlayer(Player("Kyle","Burn"))
tourney.regPlayer(Player("KG","Storm"))
tourney.regRandomPlayers(1, data)

tourney.simulate_tournament(data)

print(tourney)
'''

'''
decks = ["G Tron", "UW Control", "Titan", "Dredge", "Jund", "Burn", "Storm"]
data = MTGData()
data.regDecks(decks)
data.set_MetaShare("UW Control", 0.0697)
data.set_MetaShare("G Tron", 0.0659)
data.set_MetaShare("Titan", 0.0386)
data.set_MetaShare("Dredge", 0.0597)
data.set_MetaShare("Jund", 0.0149)
data.set_MetaShare("Burn", 0.0336)
data.set_MetaShare("Storm", 0.0149)

data.set_MatchUp("G Tron", "G Tron", 0.5)
data.set_MatchUp("G Tron", "UW Control", 0.5)
data.set_MatchUp("G Tron", "Titan", 0.405)
data.set_MatchUp("G Tron", "Dredge", 0.529)
data.set_MatchUp("G Tron", "Jund", 0.511)
data.set_MatchUp("G Tron", "Burn", 0.413)
data.set_MatchUp("G Tron", "Storm", 0.368)

data.set_MatchUp("UW Control", "UW Control", 0.5)
data.set_MatchUp("UW Control", "Titan", 0.4)
data.set_MatchUp("UW Control", "Dredge", 0.267)
data.set_MatchUp("UW Control", "Jund", 0.6)
data.set_MatchUp("UW Control", "Burn", 0.559)
data.set_MatchUp("UW Control", "Storm", 0.268)

data.set_MatchUp("Titan", "Titan", 0.5)
data.set_MatchUp("Titan", "Dredge", 0.467)
data.set_MatchUp("Titan", "Jund", 0.525)
data.set_MatchUp("Titan", "Burn", 0.679)
data.set_MatchUp("Titan", "Storm", 0.231)

data.set_MatchUp("Dredge", "Dredge", 0.5)
data.set_MatchUp("Dredge", "Jund", 0.389)
data.set_MatchUp("Dredge", "Burn", 0.652)
data.set_MatchUp("Dredge", "Storm", 0.333)

data.set_MatchUp("Jund", "Jund", 0.5)
data.set_MatchUp("Jund", "Burn", 0.588)
data.set_MatchUp("Jund", "Storm", 0.565)

data.set_MatchUp("Burn", "Burn", 0.5)
data.set_MatchUp("Burn", "Storm", 0.745)

data.set_MatchUp("Storm", "Storm", 0.5)

result = TourneyOutcomeSimulator(1000, data)

result.regPlayer(Player("Austin","G Tron"))
#result.regPlayer(Player("Stephen","UW Control"))
#result.regPlayer(Player("Rob","Titan"))
#result.regPlayer(Player("Jon","Dredge"))
#result.regPlayer(Player("Ken","Jund"))
#result.regPlayer(Player("Kyle","Burn"))
#result.regPlayer(Player("KG","Storm"))
result.regRandomPlayers(10)

result.run_simulation()

print(result)
'''

'''
data = MTGData()
data.readMatchupRates("2019-modern-winrates.txt")
data.readMetaShare("2019-modern-meta.txt")

result = TourneyOutcomeSimulator(1000, data)

result.regPlayer(Player("Austin","Tron"))
result.regRandomPlayers(12)

result.run_simulation()

print(result)
'''

data = MTGData()
data.readMatchupRates("2019-modern-winrates.txt")
data.readMetaShare("2019-modern-meta.txt")

num_sim = 10
num_events = 100
num_players = 32
result = RandomTourneyOutcomeSimulator(num_sim, num_events, num_players, data)

result.regStaticPlayer(Player("Austin","Tron"))

result.run_simulations()

print(result)

## Python Data Reader ##
## Run python scripts using Python 3.5 or newer ##

''' Reader script:
    ----------------
    Object-Oriented approach to reading and compiling competitive match results.
    
    Author:     Austin Ladshaw
    Date:       06/04/2019
    Copyright:  This software was designed and built by Austin Ladshaw. 
                Copyright (c) 2019, all rights reserved.
'''

class GPData(object):
    def __init__(self):
        self.games = {}
        self.wins = {}
        
    def __str__(self):
        string = str(self.games) + "\n\n"
        string += str(self.wins)
        return string

    def displayDecks(self):
        for deck in self.games:
            print(deck)

    def regDeck(self, deck):
        if not deck in self.games:
            self.games[deck] = {}
            self.wins[deck] = {}

    def regResult(self, deck1, deck2, games, wins):
        if not deck2 in self.games[deck1]:
            self.games[deck1][deck2] = games
            self.wins[deck1][deck2] = wins
        else:
            self.games[deck1][deck2] += games
            self.wins[deck1][deck2] += wins

    def addGPfile(self, file):
        data = open(file,"r")
        for line in data:
            strings = line.split("\t")
            strings[2] = strings[2].split("\n")
            if strings[0] == "ZZZunclassfiableDECK" or strings[0] == "ZZZunknownDECK":
                strings[0] = "Unknown"
            if strings[2][0] == "ZZZunclassfiableDECK" or strings[2][0] == "ZZZunknownDECK":
                strings[2][0] = "Unknown"
            self.regDeck(strings[0])
            deck1 = strings[0]
            deck2 = strings[2][0]
            result = 0  #Zero is a loss (i.e., deck1 lost to deck2)
            if strings[1] == "beat":
                result = 1
            self.regResult(deck1,deck2,1,result)
        data.close()

    def addMatchfile(self, file):
        data = open(file,"r")
        start = False
        type = ""
        temp_games = {}
        temp_wins = {}
        header = []
        header_set = False
        for line in data:
            strings = line.split("\t")
            size = len(strings)
            if start == True and strings[0] == '':
                start = False
            if start == False and strings[0] != '':
                type = strings[0]
                start = True
            if start == True and strings[0] == "Decks" and header_set == False:
                for i in range(1, size-1):
                    header.append(strings[i])
                header.append(strings[size-1].split("\n")[0])
                header_set = True
            elif start == True and strings[0] != "Decks" and strings[0] != "GamesWon" and strings[0] != "GamesPlayed":
                if type == "GamesWon":
                    deck1 = strings[0]
                    temp_wins[deck1] = {}
                    for i in range(0,len(header)):
                        deck2 = header[i]
                        result = strings[i+1].split("\n")[0]
                        temp_wins[deck1][deck2] = result
                if type == "GamesPlayed":
                    deck1 = strings[0]
                    temp_games[deck1] = {}
                    for i in range(0,len(header)):
                        deck2 = header[i]
                        result = strings[i+1].split("\n")[0]
                        temp_games[deck1][deck2] = result
        data.close()
        for deck1 in temp_games:
            self.regDeck(deck1)
            for deck2 in temp_games[deck1]:
                self.regResult(deck1,deck2,int(temp_games[deck1][deck2]),int(temp_wins[deck1][deck2]))

    def recordMatchResults(self, MatchRates, MatchData):
        rates = open(MatchRates, 'w')
        data = open(MatchData, 'w')

        rate_string = "WinRates\nDecks"
        data_string = "GamesWon\nDecks"
        deck2s = []
        for deck2 in self.wins:
            data_string += "\t" + deck2
            rate_string += "\t" + deck2
            deck2s.append(deck2)
        data_string += "\n"
        rate_string += "\n"
        i = 0
        for deck1 in self.wins:
            i = 0
            for deck2 in deck2s:
                if i == 0:
                    rate_string += deck1
                    data_string += deck1
                try:
                    data_string += "\t" + str(self.wins[deck1][deck2])
                    num_wins = int(self.wins[deck1][deck2])
                    num_games = int(self.games[deck1][deck2])
                    ratio = float(num_wins)/float(num_games)
                    if ratio == 1.0 and num_games == 1:
                        ratio = 0.5
                    if ratio == 0.0 and num_games == 1:
                        ratio = 0.5
                    if ratio == 1.0 and num_games == 2:
                        ratio = 0.75
                    if ratio == 0.0 and num_games == 2:
                        ratio = 0.25
                    if ratio == 1.0 and num_games == 3:
                        ratio = 0.8
                    if ratio == 0.0 and num_games == 3:
                        ratio = 0.2
                    if ratio == 1.0 and num_games > 3:
                        ratio = 0.85
                    if ratio == 0.0 and num_games > 3:
                        ratio = 0.15
                    rate_string += "\t" + str(ratio)
                except:
                    data_string += "\t" + "0"
                    rate_string += "\t" + "0.5"
                i += 1
            data_string += "\n"
            rate_string += "\n"
        data_string += "\n"
        rate_string += "\n"

        data_string += "GamesPlayed\nDecks"
        for deck2 in self.wins:
            data_string += "\t" + deck2
        data_string += "\n"
        i = 0
        for deck1 in self.wins:
            i = 0
            for deck2 in deck2s:
                if i == 0:
                    data_string += deck1
                try:
                    data_string += "\t" + str(self.games[deck1][deck2])
                except:
                    data_string += "\t" + "0"
                i += 1
            data_string += "\n"
        data_string += "\n"

        rates.write(rate_string)
        data.write(data_string)
        rates.close()
        data.close()

## Testing ##

test = GPData()
test.addGPfile("GP-Bilbao-2019-Data.txt")
test.addMatchfile("mtg-modern-data-2018.txt")
test.recordMatchResults("2019-modern-winrates.txt","2019-modern-data.txt")

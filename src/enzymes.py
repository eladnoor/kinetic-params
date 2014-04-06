#!/usr/bin/python
#
# Except for the default python packages for python 2.5, this program requires:
# sqlite3 - a python interface for SQLite
################################################################################

import os
import sys
import urllib
import types
import re
import sqlite3
import bag

################################################################################
#                               EXCEPTIONS                                     #
################################################################################

# Description of tables:
# * note that indexed columns are marked with <>
#
# kegg_compound (cid INT, name TEXT, all_names TEXT)
#     <cid>     - the ID of the compound according to KEGG
#     name      - the name of the compound (first one appearing in KEGG)
#     all_names - a single string with all the names for this compound (; separated)
#
# kegg_name_to_cid (name TEXT, cid INT)
#     <name>    - compound name
#     cid       - the ID of the compound according to KEGG
#
# kegg_reaction (rid INT, ec TEXT, all_ec TEXT, name TEXT)
#     <rid>     - the ID of the reaction according to KEGG
#     all_ec    - a single string with all the names for this compound (; separated)
#     name      - the name of the reaction
#
# kegg_rid_to_ec (rid INT, ec TEXT)
#     rid       - the ID of the reaction according to KEGG
#     ec        - the EC number (4 digits '.'-separated)
# 
# kegg_rid_to_cid (rid INT, side INT, coefficient INT, cid INT)
#     rid       - the ID of the reaction according to KEGG
#     side      - (-1) if this compound is on the left side, (+1) on the right
#     coefficient - the stoichiometric coefficient of this compound in this reaction
#     cid       - the ID of the compound according to KEGG
#
# kegg_rid_to_numsubs (rid INT, side INT, numsubs INT)
#     <rid>     - the ID of the reaction according to KEGG
#     <side>    - (-1) if this compound is on the left side, (+1) on the right
#     numsubs   - the number of substrates for this side of this reaction
#     * there is a single index for <rid, side>
#
# kegg_module (mid INT, name TEXT)
#     <mid>     - the ID of the module in KEGG
#     name      - the name of the module in KEGG
#
# kegg_mid_to_rid (mid INT, rid INT)
#     * assocaites modules to reactions
#
#
# brenda_param (field TEXT, ec TEXT, organism INT, side INT, compound TEXT, pubid INT, value REAL)
#     field     - TN (turnover number) or KM
#     ec        - EC number (4 digits)
#     organism  - the ID of the organism (the index to the table "brenda_organism")
#     side      - the direction of the reaction (-1) - left, (0) - unknown, (+1) - right
#     compound  - the name of the compound in BRENDA
#     pubid     - the ID for the publication in BRENDA
#     value     - the value Km (in units of mM) or Kcat (in units of 1/s)
#
# brenda_organism (oid INT, name TEXT)
#     <oid>     - the ID of the organism
#     name      - the name of the organism

class Common:
    @staticmethod
    def cannonic_name(compound_name):
        """
            Change the letters to lowercase and replace all spaces and dashes
            with underscores. This might help a little with standardizing
            the compound names which are nothing but standard.
        """
        s = compound_name.lower()
        s = re.sub('[^A-Z^a-z^0-9^,\+]', '', s)
        try:    
            return unicode(s)
        except UnicodeDecodeError:
            return u"?"
            
class KeggParseException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Kegg:
    def __init__(self, comm, log_file=None):
        if (log_file != None):
            self.LOG_FILE = log_file
        else:
            self.LOG_FILE = sys.stderr

        self.COMPOUND_URL = 'ftp://ftp.genome.jp/pub/kegg/ligand/compound/compound'
        self.REACTION_URL = 'ftp://ftp.genome.jp/pub/kegg/ligand/reaction/reaction'
        self.MODULE_URL = 'ftp://ftp.genome.jp/pub/kegg/pathway/module'
        self.DB = 'res/enzymes.sqlite'

        self.COMPOUND_FILE = 'data/kegg_compound.txt'
        self.REACTION_FILE = 'data/kegg_reaction.txt'
        self.MODULE_FILE = 'data/kegg_module.txt'
        
        c = comm.cursor()
        comm.commit()
        
        if (True):
            c.execute("DROP TABLE IF EXISTS kegg_compound")
            c.execute("CREATE TABLE kegg_compound (cid INT, first_name TEXT, all_names TEXT)")
            c.execute("DROP INDEX IF EXISTS cid_idx")
            c.execute("CREATE UNIQUE INDEX cid_idx ON kegg_compound (cid)")

            c.execute("DROP TABLE IF EXISTS kegg_name_to_cid")
            c.execute("CREATE TABLE kegg_name_to_cid (name TEXT, cid INT)")
            c.execute("DROP INDEX IF EXISTS compound_name_idx")
            c.execute("CREATE UNIQUE INDEX compound_name_idx ON kegg_name_to_cid (name)")

            cid2fields_map = self.parse_kegg_file(self.COMPOUND_URL, self.COMPOUND_FILE)
            self.LOG_FILE.write("Adding the compounds into kegg_compound table ... ")
    
            name_to_cid_map = {}
            for key in sorted(cid2fields_map.keys()):
                field_map = cid2fields_map[key]
                cid = int(key[1:])
                
                all_names = u'?'
                name = u'?'

                try:
                    value_name = field_map["NAME"].strip()
                except KeyError:
                    self.LOG_FILE.write("CID " + key + " doesn't have a NAME field\n")
                    value_name = '?'

                try:
                    first_name = unicode(value_name.split(';')[0])
                    all_names = unicode(value_name.replace('\t', ''))
                except UnicodeDecodeError:
                    self.LOG_FILE.write("cannot decode ASCII string: " + field_map["NAME"] + "\n")
                
                c.execute("INSERT INTO kegg_compound VALUES(?,?,?)", (cid, first_name, all_names))
                if (not (all_names == u'?')):
                    for name in all_names.split(';'):
                        cannonic_name = Common.cannonic_name(name)
                        if (cannonic_name in name_to_cid_map):
                            name_to_cid_map[cannonic_name].append(cid)
                            self.LOG_FILE.write("Cannonic name hit: " + str(cannonic_name) + " - " + str(name_to_cid_map[cannonic_name]) + "\n")
                        else:
                            c.execute("INSERT INTO kegg_name_to_cid VALUES(?,?)", (cannonic_name, cid))
                            name_to_cid_map[cannonic_name] = [cid]

            comm.commit()
            self.LOG_FILE.write(' [DONE]\n')

        if (True):
            c.execute("DROP TABLE IF EXISTS kegg_reaction")
            c.execute("CREATE TABLE kegg_reaction (rid INT, all_ec TEXT, name TEXT)")
            c.execute("DROP INDEX IF EXISTS rid_idx")
            c.execute("CREATE UNIQUE INDEX rid_idx ON kegg_reaction (rid)")

            c.execute("DROP TABLE IF EXISTS kegg_rid_to_cid")
            c.execute("CREATE TABLE kegg_rid_to_cid (rid INT, side INT, coefficient INT, cid INT)")
            c.execute("DROP TABLE IF EXISTS kegg_rid_to_ec")
            
            c.execute("CREATE TABLE kegg_rid_to_ec (rid INT, ec TEXT)")
            c.execute("DROP INDEX IF EXISTS rid_ec_idx")
            c.execute("CREATE UNIQUE INDEX rid_ec_idx ON kegg_rid_to_ec (rid, ec)")

            c.execute("DROP TABLE IF EXISTS kegg_rid_to_numsubs")
            c.execute("CREATE TABLE kegg_rid_to_numsubs (rid INT, side INT, numsubs INT)")
            c.execute("DROP INDEX IF EXISTS rid_side_idx")
            c.execute("CREATE UNIQUE INDEX rid_side_idx ON kegg_rid_to_numsubs (rid, side)")

            rid2fields_map = self.parse_kegg_file(self.REACTION_URL, self.REACTION_FILE)
            self.LOG_FILE.write("Adding the reactions into kegg_reaction table ... ")
            for key in sorted(rid2fields_map.keys()):
                field_map = rid2fields_map[key]
                rid = int(key[1:])
                
                ec_list = field_map.get("ENZYME", "-.-.-.-").split()
                all_ec = u";".join([unicode(ec) for ec in ec_list])
                name = unicode(field_map.get("NAME", "?"))
                c.execute("INSERT INTO kegg_reaction VALUES(?,?,?)", (rid, all_ec, name))
                for ec in ec_list:
                    c.execute("INSERT INTO kegg_rid_to_ec VALUES(?,?)", (rid, unicode(ec)))                

                equation_value = field_map.get("EQUATION", "<=>")
                try:
                    (left_bag, right_bag, direction) = self.parse_reaction_formula(equation_value)
                    c.execute("INSERT INTO kegg_rid_to_numsubs VALUES(?,?,?)", (rid, -1, len(left_bag)))
                    c.execute("INSERT INTO kegg_rid_to_numsubs VALUES(?,?,?)", (rid, 1, len(right_bag)))
                    for (compound, coeff) in left_bag.iteritems():
                        cid = int(compound[1:])
                        c.execute("INSERT INTO kegg_rid_to_cid VALUES(?,?,?,?)", (rid, -1, coeff, cid))
                    for (compound, coeff) in right_bag.iteritems():
                        cid = int(compound[1:])
                        c.execute("INSERT INTO kegg_rid_to_cid VALUES(?,?,?,?)", (rid, 1, coeff, cid))
                except KeggParseException:
                    pass
                except ValueError:
                    pass
                
            comm.commit()
            self.LOG_FILE.write(' [DONE]\n')

        if (True):        
            c.execute("DROP TABLE IF EXISTS kegg_module")
            c.execute("CREATE TABLE kegg_module (mid INT, name TEXT)")
            c.execute("DROP INDEX IF EXISTS mid_idx")
            c.execute("CREATE UNIQUE INDEX mid_idx ON kegg_module (mid)")
            
            c.execute("DROP TABLE IF EXISTS kegg_mid_ec_rid")
            c.execute("CREATE TABLE kegg_mid_ec_rid (mid INT, ec TEXT, rid INT)")
            c.execute("DROP TABLE IF EXISTS kegg_mid_ec_rid_temp")
            c.execute("CREATE TABLE kegg_mid_ec_rid_temp (mid INT, ec TEXT, rid INT)")
            
            mid2fields_map = self.parse_kegg_file(self.MODULE_URL, self.MODULE_FILE)
            self.LOG_FILE.write("Adding the modules into kegg_module table ... ")
            for key in sorted(mid2fields_map.keys()):
                field_map = mid2fields_map[key]
                mid = int(key[1:])
                name = unicode(field_map.get("NAME", "?"))
                c.execute("INSERT INTO kegg_module VALUES(?,?)", (mid, name))
                
                #if ("REACTION" in field_map):
                #    reactions = field_map["REACTION"]
                #    for rid in re.findall('R([0-9]+)', reactions):
                #        c.execute("INSERT INTO kegg_mid_to_rid VALUES(?,?)", (mid, int(rid)))
                if ("ORTHOLOGY" in field_map):
                    orthology = field_map["ORTHOLOGY"]
                    ec_list = []
                    rid_list = []
                    for (ec_clause, rn_clause) in re.findall('\[EC:([^\]]+)\]\s+\[RN:([^\]]+)\]', orthology):
                        ec_list = ec_clause.split()
                        rid_list = [int(rid[1:]) for rid in rn_clause.split()]
                        for ec in ec_list:
                            for rid in rid_list:
                                # now search for the RID, EC combination in kegg_reaction
                                c.execute("INSERT INTO kegg_mid_ec_rid_temp VALUES(?,?,?)", (mid, ec, rid))

            comm.commit()
            
            c.execute("INSERT INTO kegg_mid_ec_rid SELECT a.mid, a.ec, a.rid FROM kegg_mid_ec_rid_temp a, kegg_rid_to_ec b where a.ec=b.ec and a.rid=b.rid;")
            c.execute("DROP TABLE kegg_mid_ec_rid_temp")
            
            self.LOG_FILE.write(' [DONE]\n')
        c.close()
                  
    def parse_kegg_file(self, url, filename):
        if (not os.path.exists(filename)):
            self.LOG_FILE.write("Downloading from: " + url + " to " + filename + " ... ")
            urllib.urlretrieve(url, filename)
            self.LOG_FILE.write("[DONE]\n")

        self.LOG_FILE.write("Parsing file: " + filename + " ")
        kegg_file = open(filename, 'r')
        curr_field = ""
        field_map = {}
        line = kegg_file.readline()
        line_counter = 0
        entry2fields_map = {}
        while (line):
            field = line[0:12].rstrip()
            value = line[12:].strip()
    
            if (field == "///"):
                entry = field_map["ENTRY"].split()[0]
                entry2fields_map[entry] = field_map
                field_map = {}
            else:
                if (field != ""):
                    curr_field = field
                if (curr_field in field_map):
                    field_map[curr_field] = field_map[curr_field] + "\t" + value
                else:
                    field_map[curr_field] = value
    
            line = kegg_file.readline()
            line_counter += 1
            if (line_counter % 20000 == 0):
                sys.stderr.write('.')
                sys.stderr.flush()
        
        kegg_file.close()
        self.LOG_FILE.write(" [DONE]\n")
        return entry2fields_map
        
    def parse_reaction_formula_side(self, s):
        """ parse the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
            return the set of CIDs, ignore stoichiometry
        """
        compound_bag = bag.bag()
        for member in s.split("+"):
            member = member.strip()
            if (member.find(' ') > -1):
                (amount, cid) = member.split()
            else:
                amount = 1
                cid = member
            
            try:
                compound_bag.add(cid, int(amount))
            except ValueError:
                raise KeggParseException("Non-specific reaction: " + s)
        
        return compound_bag
        
      
    def parse_reaction_formula(self, formula):
        """ parse a two-sided formula such as: C00001 => C00002 + C00003 
            return the set of substrates, products and the direction of the reaction
        """
        sindex = formula.find('=')
        if (sindex == -1):
            return None
        left_side = formula[0:sindex-1].strip()
        right_side = formula[sindex+2:].strip()
        direction = formula[(sindex-1):(sindex+2)].strip() # the direction: <=, => or <=>
        
        left_bag = self.parse_reaction_formula_side(left_side)
        right_bag = self.parse_reaction_formula_side(right_side)
        
        return (left_bag, right_bag, direction)

class BrendaParseException(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

class Brenda:
    def __init__(self, comm, log_file=None):
        if (log_file != None):
            self.LOG_FILE = log_file
        else:
            self.LOG_FILE = sys.stderr
        
        self.BRENDA_FILE = 'data/brenda_fixed.txt'
        if (not os.path.exists(self.BRENDA_FILE)):
            # Fix the lines of the downloaded file (one single line per entry - no \n in the middle)
            brenda = open('data/brenda_download.txt', 'r')
            brenda_fixed = open('data/brenda_fixed.txt', 'w')
            line_number = 0
            concat_line = ""
            while (True):
                line_number += 1
                line = brenda.readline()
                if (line == ''):
                    break
                if (line[0] == '\t'):
                    concat_line += ' ' + line.strip()
                else:
                    if (concat_line != ""):
                        brenda_fixed.write(concat_line + "\n")
                    concat_line = line.strip()
            
            brenda.close()
            brenda_fixed.close()

        # Parse the fixed brenda file and add the data to the sqlite database file
        c = comm.cursor()
        c.execute("DROP TABLE IF EXISTS brenda_param")
        c.execute("CREATE TABLE brenda_param (field TEXT, ec TEXT, organism TEXT, compound TEXT, pubid INT, value REAL)")

        c.execute("DROP TABLE IF EXISTS brenda_organism")
        c.execute("CREATE TABLE brenda_organism (oid INT, name TEXT)")
        c.execute("DROP INDEX IF EXISTS oid_idx")
        c.execute("CREATE UNIQUE INDEX oid_idx ON brenda_organism (oid)")        

        comm.commit()

        self.LOG_FILE.write("Parsing the BRENDA data file ")
        enzyme_counter = 0

        brenda_file = open(self.BRENDA_FILE, 'r')
        organism_map = {}
        while (True):
            datamap = self.parse_brenda_enzyme(brenda_file)
            if (datamap == {}):
                break
            
            enzyme_counter += 1
            if (enzyme_counter % 1000 == 0):
                sys.stderr.write('.')
                sys.stderr.flush()

            if (len(datamap['ID']) != 1):
                print datamap['ID']
                raise BrendaParseException("There isn't one single ID field for enzyme #%d" % enzyme_counter)
            ec_number = datamap['ID'][0]
            
            if (len(datamap['RN']) == 0): # REACTION
                raise BrendaParseException("There isn't no RN field for enzyme #%d" % enzyme_counter)
            recommended_name = datamap['RN'][0]

            if (not 'RE' in datamap): # REACTION
                #self.LOG_FILE.write("Warning: there isn't one single RE field for enzyme #%d\n" % enzyme_counter)
                substrates = []
                products = []
            else:
                reaction = re.sub('\([^\(^\)]+\)$', '', datamap['RE'][0], count=1)
                (substrates, products) = self.parse_formula(reaction)
    
            for value in datamap.get('PR', []): # Protein
                tokens = re.split('^#(\d+)# ([A-Za-z0-9 \.]+)', value)
                if (len(tokens) != 4):
                    self.LOG_FILE.write("Warning: problem with PR line - " + value + "\n")
                    continue
                organism_map[int(tokens[1])] = unicode(tokens[2].strip())

            # example: TN	#16# 3.3 {Dihydroxyacetone}  (#16# pH 7.0, 25'C <14>) <14>
            for (field, organism_id, k, cannonic_name, pubid) in self.parse_params(datamap, ['KM', 'TN']):
                organism = organism_map[organism_id]
                c.execute("INSERT INTO brenda_param VALUES(?,?,?,?,?,?)", (field, ec_number, organism, cannonic_name, pubid, k))

        brenda_file.close()
        comm.commit()
        c.close()
        self.LOG_FILE.write("[DONE]\n")               

    def parse_brenda_enzyme(self, brenda_file):
        datamap = {}
        
        while True: # marks the end of the current enzyme
            line = brenda_file.readline().strip()
            if (line == '' or line == '///'): # End Of File
                return datamap
                
            if (len(line.split()) == 1): # there is only one token in this line, so it must be a header
                continue

            (field, value) = line.split(None, 1)
            field = field.rstrip()
            if (not field in datamap):
                datamap[field] = []
            datamap[field].append(value)

    def parse_side(self, s):
        compound_list = []
        for c in s.split(' + '):
            c = self.strip_stoichiometric_number(c)
            if (c == 'nad(p)h'):
                compound_list += ['nadh', 'nadph']
            elif (c == 'nad(p)+'):
                compound_list += ['nad+', 'nadp+']
            else:
                compound_list.append(c)
        return compound_list

    def strip_stoichiometric_number(self, s):
        count = 1
        compound = s
        tokens = s.split(' ', 1)
        if (len(tokens) == 2):
            try:
                count = int(tokens[0])
                compound = tokens[1]
            except ValueError:
                pass
        return Common.cannonic_name(compound)

    def parse_formula(self, s):
        sides = s.split(' = ', 1)
        if (len(sides) < 2):
            return ([], [])
        else:
            return (self.parse_side(sides[0]), self.parse_side(sides[1]))

    def parse_params(self, datamap, field_list):
        results = []
        for field in field_list:
            for value in datamap.get(field, []):
                if (value.find('mutant') != -1):
                    continue # skip all lines that are referring to mutants
                tokens = re.split("^#([\d,]+)#\s+([e\d\.-]+)\s+{(.*)}.+<([\d,]+)>", value)
    
                if (len(tokens) != 6):
                    self.LOG_FILE.write("Warning: problem with " + field + " line - " + value + "\n")
                    continue
    
                organism_id = int(tokens[1].split(',')[0])
    
                if (tokens[2] == '-999'):
                    continue
                k_range = tokens[2].split('-', 1)
                try:
                    k = float(k_range[-1])
                except ValueError:
                    self.LOG_FILE.write("Warning: Can't parse this value for %s - %s\n" % (field, tokens[2]))
                    continue
    
                if (tokens[3] == '' or tokens[3] == 'More'):
                    continue
                cannonic_name = Common.cannonic_name(tokens[3])
    
                pubid = int(tokens[4].split(',')[0])
                
                results.append((field, organism_id, k, cannonic_name, pubid))
        return results
        
###################################################################################################
#                                             MAIN                                                #
###################################################################################################

try:
    os.mkdir('res')
except OSError:
    pass
comm = sqlite3.connect('res/enzymes.sqlite')
KEGG = Kegg(comm)
BRENDA = Brenda(comm)

# Now Join the databases:
c = comm.cursor()
c.execute("DROP TABLE IF EXISTS merged_km_temp;")
c.execute("CREATE TABLE merged_km_temp (ec TEXT, organism TEXT, cid INT, pubid INT, value REAL);")
c.execute("INSERT INTO merged_km_temp SELECT a.ec, a.organism, b.cid, a.pubid, a.value from brenda_param a INNER JOIN kegg_name_to_cid b where a.compound=b.name and a.field='KM';");
comm.commit()

c.execute("DROP TABLE IF EXISTS merged_km;")
c.execute("CREATE TABLE merged_km (rid INT, ec TEXT, side INT, cid INT, organism TEXT, pubid INT, value REAL);")
c.execute("INSERT INTO merged_km SELECT r2e.*, r2c.side, r2c.cid, k.organism, k.pubid, k.value FROM kegg_rid_to_cid r2c, kegg_rid_to_ec r2e, merged_km_temp k where r2c.cid=k.cid and r2c.rid=r2e.rid and r2e.ec=k.ec;")

c = comm.cursor()
c.execute("DROP TABLE IF EXISTS merged_tn_temp;")
c.execute("CREATE TABLE merged_tn_temp (ec TEXT, organism TEXT, cid INT, pubid INT, value REAL);")
c.execute("INSERT INTO merged_tn_temp SELECT a.ec, a.organism, b.cid, a.pubid, a.value from brenda_param a INNER JOIN kegg_name_to_cid b where a.compound=b.name and a.field='TN';");
comm.commit()

c.execute("DROP TABLE IF EXISTS merged_tn;")
c.execute("CREATE TABLE merged_tn (rid INT, ec TEXT, side INT, cid INT, organism TEXT, pubid INT, value REAL);")
c.execute("INSERT INTO merged_tn SELECT r2e.*, r2c.side, r2c.cid, k.organism, k.pubid, k.value FROM kegg_rid_to_cid r2c, kegg_rid_to_ec r2e, merged_tn_temp k where r2c.cid=k.cid and r2c.rid=r2e.rid and r2e.ec=k.ec;")

c.close()
comm.close()

from random import uniform, randint
class Entity(object):
    def __init__(self, filename=False, alien=False, parents=False, mutant=False):
        if filename is False:
            self.genotype=zeros((NOOFRINGS,4), 'f')
            self.origin="zero"
            self.fenotype=False
            self.fitness=False
        else:
            self.genotype=load(filename)
            self.origin="zero"
            self.evaluate()
        if alien is True: self.alien()
        if parents is not False:
            (ma, pa)=parents
            self.crossover(ma, pa)
        if mutant is not False:
            self.mutant(mutant)

    def Fenotype(self):
        """Returns effective aree (the fenotype) for the given entity."""
        self.fenotype=array([array([aree[m][r]*self.genotype[:,m][r]/noofxtals[r] for r in range(NOOFRINGS)]).sum(axis=0) for m in range(3)])

    def Fitness(self):
        if self.fenotype is False: self.Fenotype()
        self.fitness=0.
        area=self.fenotype.sum(0)
        iterlist=[(i, e) for i,e in enumerate(energies) if 210.<e<511.]
        average_area=sum([area[i] for (i,e) in iterlist])/len(iterlist)
        self.fitness+=average_area**2
        for i,e in iterlist:
            #if area[i]<=350.: self.fitness-=(area[i]-1000.)**40
            #elif area[i]<=400.:
                #self.fitness-=(area[i]-1000.)**30
            #elif area[i]<=550.:
                #self.fitness+=area[i]*1.e6
            diff=abs(area[i]-area[i+1])+abs(area[i]-area[i-1])
            if diff>2.: self.fitness-=1.*(diff)**10 # condizione sulla derivata
            else: self.fitness+=1.*(diff) # condizione sulla derivata
            #elif area[i]<=600.: self.fitness-=area[i]**40
            #else: self.fitness-=area[i]**30
            self.fitness-=(area[i]-average_area)**2

    def randomgene(self, genenumber):
        """Returns a random gene"""
        probs=[uniform(0.,1.) for i in range(3)]
        # The if statements are workaround to make faster convergency,
        # by maximizing the chance that a randomgene will give a good
        # result for the fitness function
        if genenumber <26: probs[1]=0. # Rings without Cu200
        if genenumber <14: probs[0]=0. # Rings without Cu111
        probssum=sum(probs)
        gene=[int(p/probssum*noofxtals[genenumber]) for p in probs]
        gene.append(noofxtals[genenumber]-sum(gene))
        return gene

    def evaluate(self):
        self.Fenotype()
        self.Fitness()

    def alien(self):
        self.genotype=array([self.randomgene(i) for i in range(NOOFRINGS)])
        self.evaluate()
        self.origin="alien"

    def output(self):
        if self.fitness is False: self.evaluate()
        genosum=self.genotype.cumsum(axis=1)
        for i, (Z,hkl) in enumerate(matlist):
            prefix_=prefix(Z, hkl)
            save("%sarea.dat" % prefix_, array((energies, self.fenotype[i])).transpose())
            save("%s.dat" % prefix_, array((range(1,NOOFRINGS+1), genosum[:,i])).transpose())
        # Computing and saving total area
        areatot=array((energies, self.fenotype.sum(axis=0)))
        save("areatot.dat", areatot.transpose())
        save("dna.dat", self.genotype)

    def crossover(self, ma, pa):
        """mix two genotypes of the parents"""
        start= randint(14, NOOFRINGS-2)
        stop = randint(start, NOOFRINGS-1)
        for i in range(NOOFRINGS):
            if start<i<stop: self.genotype[i]=ma.genotype[i]
            else: self.genotype[i]=pa.genotype[i]
        self.evaluate()
        self.origin="offspring"

    def mutegene(self, genenumber):
        return gene+mutation

    def mutant(self, original):
        """Random mutation in only one gene: copies old string and then changes only one gene"""
        mutantgene=randint(14,NOOFRINGS-1)
        for g in range(NOOFRINGS):
            if g!=mutantgene:
                self.genotype[g]=original.genotype[g]
            else:
                gene=original.genotype[g]
                mutation=zeros(4)
                for m in range(flags['mute']):
                    # Donor selection
                    choices=[choice for choice in range(4) if (gene[choice]+mutation[choice])!=0]
                    if len(choices)>1: a = choices.pop(randint(0,len(choices)-1))
                    else: a=choices[0]
                    # Acceptor selection
                    choices=[choice for choice in range(4) if choice!=a]
                    b = choices.pop(randint(0,len(choices)-1))
                    # Calculating mutation
                    mutation[a]-=1
                    mutation[b]+=1
                self.genotype[g]=gene+mutation
                delta_fenotype=array([array([aree[material][g]*mutation[material]/noofxtals[g]]).sum(axis=0) for material in range(3)])
                self.fenotype=original.fenotype+delta_fenotype
        self.Fitness()
        self.origin="mutant"

def report(e1, e2):
    if e1.fitness<e2.fitness:
        print "EVOLUTION from a %s!!!!!!!!!!! :DDDDDD" % e2.origin
        percentual_evolution=-200.*float(e1.fitness-e2.fitness)/(e1.fitness+e2.fitness)
        print "%g\t%g\t%g\n" % (e1.fitness, e2.fitness, percentual_evolution)
    elif e1.fitness==e2.fitness:
        print "Everything still... :|"
    else:
        print "%g => %g\n" % (e1.fitness, e2.fitness)
        print "Involution from a %s T_T " % e2.origin
        exit()

def find_the_best_in(population):
    """Finds out and returns the most fitting element in a population"""
    values=[e.fitness for e in population]
    best_index=values.index(max(values))
    return population[best_index]

def evolve(population, oldbest):
    if oldbest.fitness>0.:
        threshold=oldbest.fitness*(1.-flags['gap'])
    else:
        threshold=oldbest.fitness*(1.+flags['gap'])
    survived=0
    dead=0
    offspring=0
    mutants=0
    stranger=0
    newpopulation=[]
    for entity in population:
        if entity.fitness >= threshold:
            newpopulation.append(entity)
            survived+=1
        else:
            dead+=1
    if survived==flags['population']:
        if flags['gap']<1.e-12:
            flags['gap']=0.01
            flags['population']=flags['population']*2
            if flags['population']>1.e4:
                print "\n\nGiving up with a population of %d elements" % (flags['population']/2)
                sys.exit()
            print "Catastrophe coming!!!!!!"
            survived=1
            dead=flags['population']-1
            newpopulation=[oldbest]
        else:
            print "Recalculating gap!!!!!!"
            flags['gap']=0.5*flags['gap']
        print "New gap = %g" % flags['gap']
        mutants=dead
    elif survived==1:
        mutants=dead
        stranger=dead-mutants
    else:
        mutants=randint(dead/2, dead)
        offspring=dead-mutants
        stranger=dead-offspring-mutants

    # Population changing
    for i in range(mutants):
        original=population[randint(0,survived-1)]
        newpopulation.append(Entity(mutant=original))
    for i in range(offspring):
        # Chosing parents randomly between survived population
        r1=randint(0, len(newpopulation)-2)
        r2=randint(r1, len(newpopulation)-1)
        ma, pa = population[r1], population[r2]
        newpopulation.append(Entity(parents=(ma, pa)))
    for i in range(stranger): newpopulation.append(Entity(alien=True))
    if not (survived+stranger+offspring+mutants)==flags['population']:
        print survived,stranger,offspring,mutants, survived+stranger+offspring+mutants
        exit()

    # Evaluating generation
    newbest=find_the_best_in(population)

    print "Survived: %s" % survived
    print "Strangers: %s" % stranger
    print "Offspring: %s" % offspring
    print "Mutants: %s" % mutants
    report(oldbest, newbest)
    # Returns
    return newpopulation, newbest

def init_population():
    # Empty list of genotypes
    population=[]
    # Deciding if the last best one should be included
    if flags['reload']=='ask':
        takelast=raw_input("Do you want to reload old string.dat? ")
    elif flags['reload'] is True:
        takelast='yes'
    elif flags['reload'] is False:
        takelast='no'
    if takelast.lower() in ('y', 'yes'):
        population.append(Entity("dna.dat"))
    # Including other invited guys
    [population.append(Entity(arg)) for arg in args]
    remaining_guys=flags['population'] - len(population)
    if remaining_guys<0.: flags['population']=len(population)
    else: [population.append(Entity(alien=True)) for i in range(remaining_guys)]
    return population

def dance_of_life():
    # first generation
    generation=0
    population=init_population()
    oldbest=find_the_best_in(population)
    file("generation.dat",'w').write("%s\t%s\n" % (generation, oldbest.fitness))
    while True:
        generation+=1
        myfile=file("generation.dat",'a')
        print "Generation: %s" % generation
        population, newbest = evolve(population, oldbest)
        if oldbest.fitness<newbest.fitness:
            newbest.output()
            myfile.write("%s\t%s\n" % (generation, newbest.fitness))
            myfile.close()
            oldbest=find_the_best_in(population)
        print ""

dance_of_life()

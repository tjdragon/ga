package org.tj.ga

case class NbChromos(val nb: Int) extends AnyVal
case class NbGenes(val nb: Int) extends AnyVal

abstract class GA[T, U] {
  import scala.language.implicitConversions
  implicit def nbgenes2int(obj: NbGenes): Int = obj.nb
  implicit def nbchromos2int(obj: NbChromos): Int = obj.nb

  protected val rnd = new scala.util.Random
  protected def nextInt(upper: Int) = rnd.nextInt(upper)

  type Gene = T
  type Chromo = Seq[Gene]
  type Population = Seq[Chromo]

  def ~~(): Gene // Generates a random gene
  def ~~(gene: Gene): Gene = ~~() // Mutates a given gene
  def xx(c1: Chromo, c2: Chromo): (Chromo, Chromo) = { // Crosses two chromosomes
    require(c1.length == c2.length, s"Crossing chromos of different lengths ${c1.length} != ${c2.length}")
    val idx = nextInt(c1.length)
    val rmi = c1.length - idx
    (c1.take(idx) ++ c2.takeRight(rmi), c2.take(idx) ++ c1.takeRight(rmi))
  }
  def ~~(c: Chromo): Chromo = { // Mutates a chromosome
    val idx = nextInt(c.length)
    val mc = ~~(c(idx))
    c.take(idx) ++ List(mc) ++ c.takeRight(c.length - 1 - idx)
  }
  def ~#(nbGenes: NbGenes): Chromo = // Generates a random chromosome 
    (for (i <- 0 until nbGenes) yield ~~())
  def ~#(nbChromos: NbChromos, nbGenes: NbGenes): Population = { // Generates a random population
    val nbc = if (nbChromos % 2 == 0) nbChromos.nb else 1 + nbChromos.nb // Makes sure we have an even number
    for (i <- 0 until nbc) yield ~#(nbGenes)
  }

  // The algorithm uses a 'decoder' (to decoded a chromo into the target),
  // an 'evaluator' to evaluate a solution towards its target
  // and a generic 'solve' method
  def decoder(chromo: Chromo): U
  def evaluator(solution: U, target: Option[U]): Double // The higher the number the worst, 0 is the best

  def solve(decode: Chromo => U)(evaluate: (U, Option[U]) => Double)(nbChromos: NbChromos, nbGenes: NbGenes, target: Option[U], crossOverRate: Double = 0.7, mutationRate: Double = 0.2, maxIterations: Double = 10000): (Chromo, U, Double) = {
    val nbBests = 2
    val initialPopulation = ~#(nbChromos, nbGenes)

    def ~!#(in: Population, out: Population = Seq()): Population = {
      val rndDouble = rnd.nextDouble

      if (in.isEmpty) out else {
        val c1h = in.head
        val rcs = in.tail
        val c2h = rcs.head
        if (rndDouble < mutationRate) { // Let's mutate
          val mc1 = ~~(c1h)
          val mc2 = ~~(c2h)
          ~!#(rcs.tail, List(mc1, mc2) ++ out)
        } else if (rndDouble < crossOverRate) {
          val cxs = xx(c1h, c2h)
          ~!#(rcs.tail, List(cxs._1, cxs._2) ++ out)
        } else
          ~!#(rcs.tail, List(c1h, c2h) ++ out)
      }
    } // crosses, mutates, or copies the population into a new population

    def parSort(population: Population): Population = {
      val parPopulation = population.par
      val evaledParPopulation = parPopulation.map(e => (e, evaluate(decode(e), target))).toList
      val sortedParPopulation = evaledParPopulation.sortBy(e => e._2).map(e => e._1)
      sortedParPopulation
    }

    def solve(population: Population, iter: Int = 0): Chromo = {
      // Serial version Slow (60 secs)
//      val sortedChromos = population.sortBy(w => evaluate(decode(w), target))
      
      // Parallel version Fast (5 secs)
      val sortedChromos = parSort(population)

      if (iter > maxIterations) {
        sortedChromos(0)
      } else {
        val bests = sortedChromos.slice(0, nbBests)
        val best = bests(0)
        val bestDecoded = decode(best)
        val evaled = evaluator(bestDecoded, target)
        if (iter % 1000 == 0) println(s"Current iter $iter, found $bestDecoded")
        if (evaled == 0) {
          println(s"Found a solution in $iter iterations")
          best
        } else {
          val shuffledChromos = rnd.shuffle(sortedChromos.take(sortedChromos.length - nbBests))
          val nextChromos = ~!#(shuffledChromos)
          val newPopulation = bests ++ nextChromos
          solve(newPopulation, iter + 1)
        }
      }
    } // solves the problem

    val best = solve(initialPopulation)
    val bestDecoded = decode(best)
    val bestScore = evaluate(bestDecoded, target)
    (best, bestDecoded, bestScore)
  } // def solve
}
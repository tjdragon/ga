package org.tj.ga

class FormulaFinder extends GA[Char, String] {
  private val digits = List('0', '1', '2', '3', '4', '5', '6', '7', '8', '9')
  private val operands = List('+', '-', '/', '*')
  private val domain = digits ++ operands
  private def rndElem = domain(rnd.nextInt(domain.length))
  
  override def ~~(): Gene = rndElem
  override def decoder(chromo: Chromo) = chromo.map(g => g.toString()).mkString
  override def evaluator(sol: String, tgt: Option[String]): Double = {
    try {
      val eval = MVEL.eval(sol+"+0.0")
      val tgtSol = eval.asInstanceOf[Double]
      val tgtAsD = tgt.get.toDouble
      Math.abs(tgtSol-tgtAsD)
    } catch {
      case t: Throwable => Double.MaxValue
    }
  }
}

object FormulaFinder {
  def main(args: Array[String]): Unit = {
    val sga = new FormulaFinder
    val target = "123456"
    val result = sga.solve(sga.decoder)(sga.evaluator)(NbChromos(100), NbGenes(9), Some(target))
    println(s"result: $result")
  }
}
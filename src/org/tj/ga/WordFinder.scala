package org.tj.ga

class PhraseGA extends GA[Char, String] {
  private val from = 'a'; val to = 'z'
  private def rndChar: Char = (rnd.nextInt(to-from+1)+from).toChar
  
  override def ~~(): Gene = rndChar
  override def decoder(chromo: Chromo) = chromo.map(g => g.toString).mkString
  override def evaluator(sol: String, tgt: Option[String]): Double =
    if (sol.length() != tgt.get.length()) Double.MaxValue else
      (for(i <- 0 until sol.length) yield if (sol(i).equals(tgt.get(i))) 0 else 100).sum
}

object WordFinder {
  def main(args: Array[String]): Unit = {
    val sga = new PhraseGA
    val target = "abcdefghijklmnopqrstuvwxyz"
    val targetLength = target.length
    val result = sga.solve(sga.decoder)(sga.evaluator)(NbChromos(100), NbGenes(targetLength), Some(target))
    println(s"result: $result")
  }
}
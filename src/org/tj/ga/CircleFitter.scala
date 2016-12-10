package org.tj.ga

case class Circle(val x: Int, val y: Int, val radius: Int) {
  val dxy = (x - radius, y - radius)
  val dwh = (radius * 2, radius * 2)
  val surface = Math.PI * radius.toDouble * radius.toDouble
  
  def dist(other: Circle) = Math.sqrt((x - other.x) * (x - other.x) + (y - other.y) * (y - other.y))
  def intersect(other: Circle) = (radius + other.radius) > dist(other)
  
  def draw(g2d: java.awt.Graphics2D) = g2d.drawOval(dxy._1, dxy._2, dwh._1, dwh._2)
  override def toString = s"C($x, $y, $radius)"
}

case class Area(val w: Int, val h: Int, val max: Int) {
  val maxVal = Math.max(w, h)
  val rnd = new scala.util.Random
  val circles = for(i <- 0 until max) yield Circle(ni(w), ni(h), ni((w*.2).toInt))
  val nbCircles = circles.length
  
  def ni(i: Int) = rnd.nextInt(i)
  def persist2file(best: Circle) {
    val image = new java.awt.image.BufferedImage(w, h, java.awt.image.BufferedImage.TYPE_INT_ARGB)
    val g2d = image.createGraphics
    g2d.setColor(java.awt.Color.BLACK)
    circles.foreach(c => {
      c.draw(g2d)
      g2d.drawString(c.toString, c.x, c.y)
    })
    g2d.setColor(java.awt.Color.RED)
    best.draw(g2d)
    g2d.drawString(best.toString, best.x, best.y)
    javax.imageio.ImageIO.write(image, "png", new java.io.File("area.png"))
    println("Drawn circles")
  }
}

class CircleFitter(val area: Area) extends GA[Int, Circle] {
  private def rndVal = rnd.nextInt(area.maxVal)
  private val maxRadius = Math.min(area.h/2, area.w/2)
  private val maxSurface = Math.PI * maxRadius * maxRadius
  
   override def ~~(): Gene = rndVal
   override def decoder(chromo: Chromo) = Circle(chromo(0), chromo(1), chromo(2))
   override def evaluator(sol: Circle, tgt: Option[Circle]): Double = {
    // bigger the surface, the better
    val surfaceScore = Math.abs(1.0 - sol.surface / maxSurface)
    // fewer number of intersections, the better
    val nbIntersects = (for(c <- area.circles) yield if (sol.intersect(c)) 1.0 else 0.0).sum
    val interesectScore = nbIntersects / area.nbCircles
    // the solution must be strongly inside
    val insideScore = if ((sol.x-sol.radius < 0) || (sol.x+sol.radius > area.w) ||
                          (sol.y-sol.radius < 0) || (sol.y+sol.radius > area.h)) 10.0 else 0.0
    surfaceScore + interesectScore * 3.0 + insideScore
  }
}

object CircleFitter {
  def main(args: Array[String]): Unit = {
    val area = Area(800, 600, 25)
    val circleFitter = new CircleFitter(area)
    val best = circleFitter.solve(circleFitter.decoder)(circleFitter.evaluator)(NbChromos(1000), NbGenes(3), None, .7, .2, 1000)
    area.persist2file(best._2)
  }
}
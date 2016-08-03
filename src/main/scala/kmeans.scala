//import scala.collection.mutable.Seq
import scala.collection.parallel.ParSeq
import scala.io.Source
import scala.collection.mutable.ArrayBuffer

object Main {
  def main(args: Array[String]){
    //test
    if (args.length > 0){
      val q = input(args(0))
      if(args.length > 2){
        val i = args(2).toInt
          val qmat = q.map{_.take(i)}
            (2 to args(1).toInt).foreach{j =>
              bench_iter(j, q.size, q(1).size, qmat, true)
            }
      }else{
        (2 to args(1).toInt).foreach{i =>
          bench_iter(i, q.size, q(1).size, q, true)
        }
      }
    }else{
      benchmark
    }
  }

  def test{
    val rand = new scala.util.Random(1)
    //val query : Seq[Seq[Double]] = Seq.fill(10)(Seq.fill(8)(rand.nextDouble))
    val query: Seq[Seq[Double]] =
      Seq(
        Seq(4,5),
        Seq(5,5),
        Seq(4,4),
        Seq(0,5),
        Seq(3,2),
        Seq(0,4)
      )
    val clus = 2
    val kmeans = new Kmeans(query, clus)
    kmeans.repeatCrusterCalc
    println(kmeans.recruster)

    val hamerly = new Hamerly(query, clus)
    while(hamerly.eval > 0){}
    println(hamerly.a)

    println(kmeans.recruster sameElements hamerly.a)
  }

  def input(file: String):Seq[Seq[Double]] ={
    var query: Seq[Seq[Double]] = Seq.empty
    var flag: Boolean = false
    for(line <- Source.fromFile(file).getLines()){
      if(!line.startsWith("!") && !line.startsWith("\"ID_REF") && flag){
        query = query :+ line.split('\t').drop(1).map{_.toDouble}.toList
      }
      if(line.contains("!series")){
        flag = true
      }
    }
    query
  }

  def arrayCount(i:Int, array:Seq[Int])={
    (0 until i).map{num => (num, array.filter(_==num).size)}.toMap
  }

  def bench_iter(i:Int, datum:Int, dim:Int, query:Seq[Seq[Double]], opt:Boolean){
    println(i,"crusters", datum,"data", dim, "dimentions")
    printExecutionTime{
      val kmeans = new Kmeans(query, i)
      val kmsIter = kmeans.repeatCrusterCalc
      println(kmsIter, "times iter in Lloyd.", kmsIter * dim, "times point change.")
      if(opt){
        //println(kmeans.recruster)
        println(arrayCount(i, kmeans.resultCruster))
        println(kmeans.represent)
      }
    }
    printExecutionTime{
      val hamerly = new Hamerly(query, i)
      var it = 0
      var k = 0
      var change = 0
      do{
        k = hamerly.eval
        change += k
        it += 1
      }while(k>0)
        println(it, "times iter in Hamerly.", change, "times point change.")
      if(opt){
        //println(hamerly.a)
        println(arrayCount(i, hamerly.a))
      }
    }
  }

  def benchmark{
    val cruster : Seq[Int] = (1 to 3).map{scala.math.pow(10,_).asInstanceOf[Int]}
    val rand = new scala.util.Random(1)
    val data = (3 to 5).map{scala.math.pow(10, _).asInstanceOf[Int]}
    val dims = Seq(2, 10, 100)
    data.foreach{ datum =>
      dims.foreach{ dim =>
        val query : Seq[Seq[Double]] = Seq.fill(datum)(Seq.fill(dim)(rand.nextDouble))
        cruster.foreach{i =>
          bench_iter(i, datum, dim, query, false)
        }
      }
    }
  }

  def printExecutionTime(proc: => Unit) = {
    val start = System.currentTimeMillis
    proc
    println((System.currentTimeMillis - start) + "msec")
  }
}

class Hamerly(p: Seq[Seq[Double]], clusterCount: Int) extends Kmeans(p, clusterCount){
  var u = new Array[Double](p.size)
  var l = new Array[Double](p.size)
  var q = new Array[Int](clusterCount)
  var c = scala.collection.mutable.Seq.fill(clusterCount)(Seq.fill(p(0).size)(0.0))
  var a = ArrayBuffer.empty[Int]

  (0 until p.size).foreach{ i =>
    val newa = (0 until represent.size).par.map{ j => euclid_distance(p(i), represent(j))}.zipWithIndex.minBy(_._1)._2
    q(newa) = q(newa) + 1
    c(newa) = c(newa).zip(p(i)).map{e => e._1 + e._2}
    a += newa
  }
  a.toArray
  
  var cc = c.zipWithIndex.par.map{cj => cj._1.map{div(_, q(cj._2))}}
  (0 until p.size).par.foreach{i =>
    u(i) = euclid_distance(p(i), cc(a(i)))
    l(i) = cc.zipWithIndex.filter(_._2 != a(i) ).map{ci => euclid_distance(p(i), ci._1)}.min
  }

  def eval : Int = {
    val s = cc.zipWithIndex.par.map{cj => cc.zipWithIndex.filter(_._2!=cj._2).par.map{cjp => euclid_distance(cjp._1, cj._1)}.min}
    (0 until p.size).foreach{i =>
      val m = s(a(i))/2 max l(i)
      if(u(i) >= m){
        u(i) = euclid_distance(p(i), cc(a(i)))
        if(u(i) >= m){
          val aorg = a(i)
          a(i) = cc.par.map{cci => euclid_distance(p(i), cci)}.zipWithIndex.minBy(_._1)._2
          if (aorg != a(i)){
            q(aorg) -= 1
            q(a(i)) += 1
            c(aorg) = c(aorg).zip(p(i)).map{e => e._1 - e._2}
            c(a(i)) = c(a(i)).zip(p(i)).map{e => e._1 + e._2}
            u(i) = euclid_distance(p(i), cc(a(i)))
            l(i) = cc.zipWithIndex.par.filter(_._2 != a(i)).par.map{ci => euclid_distance(p(i), ci._1)}.min
          }
        }
      }
    }
    val ccnew = c.zipWithIndex.par.map{cj => cj._1.map{div(_, q(cj._2))}}
    val pd = (0 until c.size).map{i => euclid_distance(cc(i), ccnew(i))}
    val repeatCounter = pd.par.filter(_!=0).size

    cc = ccnew

    val r = pd.zipWithIndex.maxBy(_._1)._2
    (0 until p.size).par.foreach{i =>
      u(i) += pd(a(i))
      l(i) -= pd(r)
    }
    return repeatCounter
  }

  def div(x: Double, y:Double) = if ( y != 0) x/y else 0

  def euclid_distance(p1:Seq[Double], p2:Seq[Double]): Double = {
    Math.sqrt((0 until p1.size).map{ i => (p1(i) - p2(i)) * (p1(i) - p2(i)) }.sum)
  }
}


class Kmeans(p: Seq[Seq[Double]], k: Int){
  var represent = (0 until k).par.map{ j => p(j) }
  val nil = represent(0).map{ i => 0.0 }

  def recruster : ParSeq[Int] = {
    p.par.map({ point =>
      represent.map{ rep =>
        (0 until rep.size).map{ i =>
          (rep(i) - point(i)) * (rep(i) - point(i))
        }.sum
      }.zipWithIndex.minBy(_._1)._2
    })
  }

  def resultCruster : Seq[Int] = {
    p.map({ point =>
      represent.map{ rep =>
        (0 until rep.size).map{ i =>
           (rep(i) - point(i)) * (rep(i) - point(i))
        }.sum
      }.zipWithIndex.minBy(_._1)._2
    })
  }

  def recalc(cluster: ParSeq[Int]) = {
    represent = represent.zipWithIndex.par.map({ case (item, index) =>
      cluster.zip(p).filter(_._1 == index).par.map(_._2).foldRight(nil)((x:Seq[Double], y:Seq[Double]) =>
        (0 until x.size).map{ i =>
          (x(i) + y(i))
        }).map{cp => cp / cluster.filter(_==index).size}
    })
  }

  def diff(old: ParSeq[Seq[Double]]) = {
    old.zip(represent).par.map{ oldrep =>
      (0 until p(0).size).par.map{ i => (oldrep._1(i) - oldrep._2(i)) * (oldrep._1(i) - oldrep._2(i)) }.sum
    }.sum
  }

  def repeatCrusterCalc = {
    var difference = 1.0
    var index = 0
    while(difference > 0){
      index += 1
      val old = represent
      recalc(recruster)
      difference = diff(old)
      //println(difference)
    }
    index
  }
}

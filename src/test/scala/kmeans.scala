
import org.scalacheck.Properties
import org.scalacheck.Prop.{forAll,forAllNoShrink}
import org.scalacheck.Gen.choose

object KMeansSpecification extends Properties("KMeans"){

  property("repeatCrusterCalc") = forAll(choose(2,1000), choose(2,100), choose(2,20), choose(2,20)){
    (a: Int, b: Int, c: Int, d: Int) => {
      val rand = new scala.util.Random(a)
      val query : Seq[Seq[Double]] = Seq.fill(d + b)(Seq.fill(c)(rand.nextDouble))

      val kmeans = new Kmeans(query, d)
      kmeans.repeatCrusterCalc

      val hamerly = new Hamerly(query, d)
      while(hamerly.eval>0){}

      kmeans.recruster sameElements hamerly.a
    }
  }

}

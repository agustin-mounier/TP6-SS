package IntegrationMethod;

import java.util.HashSet;
import java.util.Set;

/**
 * Created by amounier on 5/13/17.
 */
public class PedestrianForce {

    Particle p;
    private static final double kn = 1.2E5;
    private static final double kt = 2.4E5;

    private double A = 2000;
    private double B = 0.08;

    private Set<Particle> neighbours;

    private static final int W = 20;
    private static final double D = 1.2;
    private static final double MARGIN = 0.1;
    private static final int BASE_DISTANCE = 6;

//    private static final double DESIRE_SPEED = 6;


    public PedestrianForce(Particle p, Set<Particle> neighbours) {
        super();
        this.p = p;
        this.neighbours = neighbours;
    }

    public double getXForce() {
        double fXSum = 0;
        for (Particle neighbour : neighbours) {
            fXSum += getFN(p, neighbour) * getENX(p, neighbour) + getFT(p, neighbour) * (-(getENY(p, neighbour)));
            //Social force
            fXSum += (-1) * A * Math.exp((-getEpsilon(p, neighbour)) / B) * getENX(p, neighbour);
        }

        double x = getDesireXDistance();
        double y = getDesireYDistance();

        fXSum += p.mass * (((p.getDesiredSpeed() * getENX(p, x, y)) - p.velX)) / 0.5;

        return fXSum;
    }

    public double getYForce() {
        double fYSum = 0;
        for (Particle neighbour : neighbours) {
            fYSum += getFN(p, neighbour) * getENY(p, neighbour) + getFT(p, neighbour) * getENX(p, neighbour);
            //Social force
            fYSum += (-1) * A * Math.exp((-getEpsilon(p, neighbour)) / B) * getENY(p, neighbour);
        }
        double x = getDesireXDistance();
        double y = getDesireYDistance();

        fYSum += p.mass * (((p.getDesiredSpeed() * getENY(p, x, y)) - p.velY)) / 0.5;

        return fYSum;
    }

    private double getENY(Particle p1, Particle p2) {
        return (p2.y - p1.y) / getDistance(p2.x, p2.y, p1.x, p1.y);
    }

    private double getENY(Particle p1, double x, double y) {
        return (y - p1.y) / getDistance(x, y, p1.x, p1.y);
    }

    private double getENX(Particle p1, Particle p2) {
        return (p2.x - p1.x) / getDistance(p2.x, p2.y, p1.x, p1.y);
    }

    private double getENX(Particle p1, double x, double y) {
        return (x - p1.x) / getDistance(x, y, p1.x, p1.y);
    }

    private double getFN(Particle p1, Particle p2) {
        return getEpsilon(p1, p2) * -kn;
    }

    private double getFT(Particle p1, Particle p2) {
        return -kt * getEpsilon(p1, p2) * (((p1.velX - p2.velX) * (-getENY(p1, p2)))
                + ((p1.velY - p2.velY) * (getENX(p1, p2))));
    }

    private double getEpsilon(Particle p1, Particle p2) {
        double ep = p1.radius + p2.radius - (getDistance(p1.x, p1.y, p2.x, p2.y));
        if (ep > p1.radius || ep > p2.radius) {
            ep = Math.min(p1.radius, p2.radius) / 2.0;
        }
        return ep;
    }

    private double getDistance(double x0, double y0, double x1, double y1) {
        return Math.sqrt(Math.pow(x0 - x1, 2) + Math.pow(y0 - y1, 2));
    }

    public void setNeighbours(Set<Particle> neighbours) {
        Set<Particle> filteredNeighbours = new HashSet<>();
        for (Particle neighbour : neighbours) {
            if (getDistance(p.getX(), p.getY(), neighbour.getX(), neighbour.getY()) <= p.getRadius() + neighbour.getRadius()) {
                filteredNeighbours.add(neighbour);
            }
        }
        this.neighbours = filteredNeighbours;
    }

    public double getDesireXDistance() {
        if (p.x - p.getRadius() <= W / 2 - D / 2 + MARGIN)
            return W / 2 - D / 2 + p.getRadius() + MARGIN;
        if (p.x + p.getRadius() >= W / 2 + D / 2 - MARGIN)
            return W / 2 + D / 2 - p.getRadius() - MARGIN;
        return p.x;
    }

    public double getDesireYDistance() {
        if (p.y < BASE_DISTANCE)
            return 0;
        return BASE_DISTANCE;
    }



    public Particle getParticle() {
        return p;
    }

    public void reset() {
        neighbours = null;
    }
}

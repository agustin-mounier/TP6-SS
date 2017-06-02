package SiloSimulation;

import CellIndexMethod.CellIndexMethod;
import IntegrationMethod.Beeman;
import IntegrationMethod.PedestrianForce;
import IntegrationMethod.Particle;
import javafx.util.Pair;

import java.io.IOException;
import java.io.PrintWriter;
import java.util.*;
import java.util.function.Consumer;

/**
 * Created by amounier on 5/13/17.
 */
public class SiloGranularSimulation {

    double L;
    double W;
    double D;

    int N;

    static final double MASS = 80;
    static final double MAX_TRIES = 100;

    static final double BELOW_SPACE = 6;

    List<Particle> particles;
    // Cant particles that went out vs time it went out
    private static Map<Integer, Double> flowint;

    int id = 1;
    double tf;
    double dt = 0.0001;
    int fpsModule = (int) (1 / dt) / 60;

    int secondsModule = (int) (1 / dt);

    Map<Particle, Beeman> integrationMethod;

    int particlesFlowCount;
    int particlesFlowCountGlobal;


    public SiloGranularSimulation(double L, double W, double D, int maxParticles, double tf) {
        this.L = L;
        this.W = W;
        this.D = D;
        this.tf = tf;
        integrationMethod = new HashMap<>();
        flowint = new HashMap<>();
        locateParticles(maxParticles);
    }

    public void run() {
        double currentTime = 0;

        CellIndexMethod CIM = new CellIndexMethod(L + BELOW_SPACE, 1, particles, false);
        int iteration = 0;
        while (N > 0) {
            if (iteration % fpsModule == 0) {
                printState(iteration / fpsModule, particles);
            }
            CIM.reloadMatrix(particles);

            Map<Particle, Set<Particle>> neighbors = CIM.findNeighbors(particles);

            List<Particle> nextGen = new ArrayList<>();

            addWallParticles(neighbors);

            for (Particle particle : particles) {
                Beeman beeman = integrationMethod.get(particle);
                beeman.getForce().setNeighbours(neighbors.get(particle));
                Particle p = beeman.moveParticle();

                double particleTotalForce = Math.abs(beeman.getForce().getXForce()) + Math.abs(beeman.getForce().getYForce());

                p.setForce(particleTotalForce);

                if (p.getY() < BELOW_SPACE && !p.isOut()) {
                    p.setOut(true);

                    particlesFlowCount++;
                    particlesFlowCountGlobal++;
                    flowint.put(particlesFlowCountGlobal, currentTime);
                }
                if (p.getY() >= BELOW_SPACE / 2) {
                    nextGen.add(p);
                }

            }
            if (iteration % secondsModule == 0) {
                particlesFlowCount = 0;
            }
            this.particles = nextGen;
            this.N = nextGen.size();
            currentTime += dt;
            iteration++;
        }
//        System.out.println("Global Flow: " + particlesFlowCountGlobal / 4.0 + " particles/s");


//        System.out.println("Global Flow: " + particlesFlowCountGlobal / 5.0 + " particles/s");

    }

    private void locateParticles(int size) {
        this.particles = new ArrayList<>();

        boolean flag = true;

        double Dmin = 0.5;
        double Dmax = 0.58;

        double desireSpeedMin = 0.8;
        double desireSpeedMax = 6;

        double fixedDesireSpeed = 2;

        while (flag && particles.size() < size) {
            double diameter = Math.random() * (Dmax - Dmin) + (Dmin);
//            double desireSpeed = Math.random() * (desireSpeedMax - desireSpeedMin) + (desireSpeedMin);
            double desireSpeed = fixedDesireSpeed;
            double r = diameter / 2;
            double x, y;

            int tries = 0;
            do {
                x = Math.random() * (this.W - 2 * r) + r;
                y = BELOW_SPACE + Math.random() * (this.L - 2 * r) + r;
                tries++;
                if (tries == MAX_TRIES) {
                    flag = false;
                    break;
                }
            } while (isOverlaped(x, y, diameter / 2));
            if (flag) {
                double mass = 60 + ((diameter - 0.5) * 380);
                Particle p = new Particle(id++, x, y, 0, 0, diameter / 2, mass, desireSpeed);
                PedestrianForce gf = new PedestrianForce(p, null);
                Beeman beeman = new Beeman(gf, dt);
                integrationMethod.put(p, beeman);
                this.particles.add(p);
            }
        }
        this.N = this.particles.size();
    }

    private void relocateParticle(Particle p) {
        double x, y;
        int tries = 0;
        do {
            x = Math.random() * (this.W - 2 * p.getRadius()) + p.getRadius();
            y = BELOW_SPACE + L - (2 * p.getRadius()) + ((tries / 10) * 2 * p.getRadius());
            tries++;
        } while (isOverlaped(x, y, p.getRadius()));
        p.setX(x);
        p.setY(y);
        p.setVelX(0);
        p.setVelY(0);
        p.setOut(false);
        PedestrianForce gf = new PedestrianForce(p, null);
        Beeman beeman = new Beeman(gf, dt);
        integrationMethod.put(p, beeman);
    }

    private void removeParticle(Particle p) {
        particles.remove(p);
        integrationMethod.remove(p);
    }

    private boolean isOverlaped(double x, double y, double r) {
        for (Particle p : particles) {
            if (getDistance(p.getX(), p.getY(), x, y) < (p.getRadius() + r))
                return true;
        }
        return false;
    }

    private double getDistance(double x0, double y0, double x1, double y1) {
        return Math.sqrt(Math.pow(x0 - x1, 2) + Math.pow(y0 - y1, 2));
    }

    public void addWallParticles(Map<Particle, Set<Particle>> neighbors) {
        for (Particle p : neighbors.keySet()) {
            List<Particle> wallParticles = getWallParticle(p);
            neighbors.get(p).addAll(wallParticles);
        }
    }

    public List<Particle> getWallParticle(Particle particle) {
        List<Particle> ret = new ArrayList<>();

        Particle p = null;
        // X - right
        p = new Particle(N + 1, -particle.getRadius(), particle.getY(), 0, 0, particle.getRadius(), particle.getMass(), 0);
        p.setWall(true);
        ret.add(p);
        // X - left
        p = new Particle(N + 2, W + particle.getRadius(), particle.getY(), 0, 0, particle.getRadius(), particle.getMass(), 0);
        p.setWall(true);
        ret.add(p);

        // Y - BOTTOM
        if (particle.getX() <= (W / 2 - D / 2) || particle.getX() >= (W / 2 + D / 2)) {
            p = new Particle(N + 3, particle.getX(), BELOW_SPACE - particle.getRadius(), 0, 0, particle.getRadius(), particle.getMass(), 0);
            p.setWall(true);
            ret.add(p);
        } else {
            if (Math.abs(particle.getY() - BELOW_SPACE) < particle.getRadius()) {
                if (particle.getX() - particle.getRadius() <= (W / 2 - D / 2)
                        && getDistance(particle.getX(), particle.getY(), (W / 2 - D / 2), BELOW_SPACE) < particle.getRadius()) {
                    p = new Particle(N + 4, (W / 2 - D / 2), BELOW_SPACE, 0, 0, particle.getRadius(), particle.getMass(), 0);
                    p.setWall(true);
                    ret.add(p);
                } else if (particle.getX() + particle.getRadius() >= (W / 2 + D / 2)
                        && getDistance(particle.getX(), particle.getY(), (W / 2 + D / 2), BELOW_SPACE) < particle.getRadius()) {
                    p = new Particle(N + 5, (W / 2 + D / 2), BELOW_SPACE, 0, 0, particle.getRadius(), particle.getMass(), 0);
                    p.setWall(true);
                    ret.add(p);
                }
            }
        }


        return ret;
    }

    public void printState(int iteration, List<Particle> particles) {

        System.out.println(N);
        System.out.println(iteration);
        for (Particle p : particles) {
            System.out.println(p);
        }
    }

    private static void printFlowIntoFile(int n) {
        try {
            PrintWriter writer = new PrintWriter("flow-" + n + ".tsv", "UTF-8");
            for (Integer i : flowint.keySet()) {
                writer.println(String.format(Locale.FRENCH, "%.3f", flowint.get(i)) + "\t" + i);
            }
            writer.close();
        } catch (IOException e) {
            // do something
        }
    }

    public static void main(String[] args) {

        int L = 20;
        int W = 20;
        double D = 1.2;
        int PARTICLES = 150;
        int SIMULATION_TIME = 15;

        SiloGranularSimulation SGM = new SiloGranularSimulation(L, W, D, PARTICLES, SIMULATION_TIME);

        SGM.run();
        printFlowIntoFile(PARTICLES);
//        System.out.println("max force " + maxForce);

    }
}
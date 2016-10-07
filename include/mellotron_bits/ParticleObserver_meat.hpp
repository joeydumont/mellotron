template <class FieldModel>
ParticleObserver<FieldModel>::ParticleObserver(Particle<FieldModel> & my_particle)
: particle(my_particle)
{}

template <class FieldModel>
ParticleObserver<FieldModel>::operator() (const arma::colvec::fixed<8> &x,
                                                double                  t)
{
  // Resize the matrices.
  const int n_rows = position.n_rows;
  position.resize(n_rows+1,3);
  momentum.resize(n_rows+1,3);
  electric_field.resize(n_rows+1,3);
  magnetic_field.resize(n_rows+1,3);

  // Push the data to the appropriate containers.
  position.row(n_rows-1) = x.subvec(1,3);
  momentum.row(n_rows-1) = x.subvec(5,7);
  electric_field.row(n_rows-1) = particle.GetElectricField();
  magnetic_field.row(n_rows-1) = particle.GetMagneticField();

  gamma.push(x(4));
  times.push(x[0]);
  chi.push(particle.GetChi());
}

template <class FieldModel>
ParticleObserver<FieldModel>::OutputData()
{
  
}
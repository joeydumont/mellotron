#ifndef PARTICLE_OBSERVER_MEAT_HPP
#define PARTICLE_OBSERVER_MEAT_HPP

namespace mellotron {

template <class FieldModel>
inline
ParticleObserver<FieldModel>::ParticleObserver(Particle<FieldModel> & my_particle, const uint my_init_size)
: particle(my_particle)
, init_size(my_init_size)
, step_counter(0)
{
  position.set_size(3,init_size);
  momentum.set_size(3,init_size);
  electric_field.set_size(3,init_size);
  magnetic_field.set_size(3,init_size);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::operator() (const arma::colvec::fixed<8> &x,
                                                double                  t)
{
  // Resize the matrices.
  const int n_cols = step_counter;step_counter++;

  if (n_cols >= init_size)
  {
    position.resize(3,n_cols+1);
    momentum.resize(3,n_cols+1);
    electric_field.resize(3,n_cols+1);
    magnetic_field.resize(3,n_cols+1);
  }

  // Push the data to the appropriate containers.
  position.col(n_cols) = x.subvec(1,3);
  momentum.col(n_cols) = x.subvec(5,7);
  gamma.push_back(std::sqrt(1.0+std::pow(arma::norm(momentum.col(n_cols),2)/particle.GetMass(),2)));
  times.push_back(x[0]);

  // Compute the electromagnetic field and store it.
  particle.ComputeFieldTensor(x[0],x[1],x[2],x[3]);
  electric_field.col(n_cols) = particle.GetElectricField();
  magnetic_field.col(n_cols) = particle.GetMagneticField();

  // Compute chi.
  double pdotE         = arma::dot(momentum.col(n_cols),electric_field.col(n_cols));
  arma::colvec lorentz = gamma[n_cols]*electric_field.col(n_cols) + arma::cross(momentum.col(n_cols),magnetic_field.col(n_cols));
  double chi_prefac    = constants::physics::hbar*particle.GetUnitSystem().omega_0_SI/(particle.GetMass()*constants::physics::electron_mass*std::pow(constants::physics::c,2));
  double chi_sq        = std::pow(chi_prefac,2)*std::pow(particle.GetMass(),-4)*(std::pow(arma::norm(lorentz,2),2)-std::pow(pdotE,2));
  chi.push_back(std::sqrt(chi_sq));
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::OutputData()
{

  // We create a hash of the initial conditions to be used
  // for the HDF5 filename and main group.
  std::size_t hash_seed = 0;
  for (uint i=0; i<3; i++)
  {
    boost::hash_combine(hash_seed, position(i,0));
    boost::hash_combine(hash_seed, momentum(i,0));
  }

  const std::string particle_label = std::to_string(hash_seed)+std::string(".hdf5");

  // Creation of the HDF5 file and main group.
  hid_t file_id  = H5Fcreate(particle_label.c_str(),
                             H5F_ACC_TRUNC,
                             H5P_DEFAULT,
                             H5P_DEFAULT);

  hid_t group_id = H5Gcreate(file_id,
                             particle_label.c_str(),
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);

  WriteAllData(group_id);

  // We close the HDF5 objects.
  H5Gclose(group_id);
  H5Fclose(file_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteAllData(hid_t group_id)
{
  WriteTimes(group_id);
  WriteGamma(group_id);
  WriteChi(group_id);
  WriteStateVector(group_id);
  WriteElectromagneticField(group_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteTimes(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Creation of the temporal dataset.
  hsize_t chunk_size_1d = 5;
  hsize_t size_dataspace = times.size();
  SetHDF5Properties(dataspace_id,plist_id, 1, &size_dataspace, &chunk_size_1d, 5u);

  // Create the group and output the data.
  dataset_id   = H5Dcreate(group_id,
                           "times",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT,
                           plist_id,
                           H5P_DEFAULT);

  status       = H5Dwrite(dataset_id,
                          H5T_NATIVE_DOUBLE,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          times.data());

  H5Dclose(dataset_id);
  H5Pclose(plist_id);
  H5Sclose(dataspace_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteGamma(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Creation of the gamma dataset.
  hsize_t chunk_size_1d = 5;
  hsize_t size_dataspace = gamma.size();
  SetHDF5Properties(dataspace_id, plist_id, 1, &size_dataspace, &chunk_size_1d, 5u);

  // Creation of the gamma dataset.
  dataset_id   = H5Dcreate(group_id,
                           "gamma",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT,
                           plist_id,
                           H5P_DEFAULT);

  status       = H5Dwrite(dataset_id,
                          H5T_NATIVE_DOUBLE,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          gamma.data());

  H5Dclose(dataset_id);
  H5Pclose(plist_id);
  H5Sclose(dataspace_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteChi(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Creation of the gamma dataset.
  hsize_t chunk_size_1d = 5;
  hsize_t size_dataspace = chi.size();
  SetHDF5Properties(dataspace_id, plist_id, 1, &size_dataspace, &chunk_size_1d, 5u);

  // Creation of the chi dataset.
  dataset_id   = H5Dcreate(group_id,
                           "chi",
                           H5T_NATIVE_DOUBLE,
                           dataspace_id,
                           H5P_DEFAULT,
                           plist_id,
                           H5P_DEFAULT);

  status       = H5Dwrite(dataset_id,
                          H5T_NATIVE_DOUBLE,
                          H5S_ALL,
                          H5S_ALL,
                          H5P_DEFAULT,
                          chi.data());

  H5Pclose(plist_id);
  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteStateVector(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Creation of the position dataset.
  hsize_t size_dataspace[2] = {position.n_cols, position.n_rows};
  hsize_t chunk_size[2]     = {50,3};
  SetHDF5Properties(dataspace_id, plist_id, 2, size_dataspace, chunk_size, 5u);

  // Creation and output.
  dataset_id                 = H5Dcreate(group_id,
                                         "position",
                                         H5T_NATIVE_DOUBLE,
                                         dataspace_id,
                                         H5P_DEFAULT,
                                         plist_id,
                                         H5P_DEFAULT);

  status                     = H5Dwrite(dataset_id,
                                        H5T_NATIVE_DOUBLE,
                                        H5S_ALL,
                                        H5S_ALL,
                                        H5P_DEFAULT,
                                        position.memptr());

  H5Dclose(dataset_id);

  // Creation of the momentum dataset.
  dataset_id                 = H5Dcreate(group_id,
                                         "momentum",
                                         H5T_NATIVE_DOUBLE,
                                         dataspace_id,
                                         H5P_DEFAULT,
                                         plist_id,
                                         H5P_DEFAULT);

  status                     = H5Dwrite(dataset_id,
                                        H5T_NATIVE_DOUBLE,
                                        H5S_ALL,
                                        H5S_ALL,
                                        H5P_DEFAULT,
                                        momentum.memptr());

  H5Dclose(dataset_id);
  H5Pclose(plist_id);
  H5Sclose(dataspace_id);
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::WriteElectromagneticField(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  hsize_t size_dataspace[2] = {electric_field.n_cols, electric_field.n_rows};
  hsize_t chunk_size[2]     = {50,3};
  SetHDF5Properties(dataspace_id, plist_id, 2, size_dataspace, chunk_size, 5);

  // Creation of the electric field dataset.
  dataset_id                 = H5Dcreate(group_id,
                                         "electric_field",
                                         H5T_NATIVE_DOUBLE,
                                         dataspace_id,
                                         H5P_DEFAULT,
                                         plist_id,
                                         H5P_DEFAULT);

  status                     = H5Dwrite(dataset_id,
                                        H5T_NATIVE_DOUBLE,
                                        H5S_ALL,
                                        H5S_ALL,
                                        H5P_DEFAULT,
                                        electric_field.memptr());

  H5Dclose(dataset_id);

  // Creation of the magnetic field dataset.
  dataset_id                 = H5Dcreate(group_id,
                                         "magnetic_field",
                                         H5T_NATIVE_DOUBLE,
                                         dataspace_id,
                                         H5P_DEFAULT,
                                         plist_id,
                                         H5P_DEFAULT);

  status                     = H5Dwrite(dataset_id,
                                        H5T_NATIVE_DOUBLE,
                                        H5S_ALL,
                                        H5S_ALL,
                                        H5P_DEFAULT,
                                        magnetic_field.memptr());

  H5Dclose(dataset_id);
  H5Pclose(plist_id);
  H5Sclose(dataspace_id);

}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::SetHDF5Properties(       hid_t    &  dataspace_id,
                                                       hid_t    &  plist_id,
                                                const  int         dim,
                                                const  hsize_t  *  size,
                                                const  hsize_t  *  chunk_size,
                                                const  uint        compression_level)
{
  // We set a simple dataspace id.
  dataspace_id = H5Screate_simple(dim, size, NULL);

  // We create a property list id, then enable chunking and compression.
  plist_id     = H5Pcreate(H5P_DATASET_CREATE);
  H5Pset_chunk(plist_id,dim,chunk_size);
  H5Pset_deflate(plist_id,compression_level);
}

template <class FieldModel>
inline
ParticleObserverLienardWiechert<FieldModel>::ParticleObserverLienardWiechert(       Particle<FieldModel>  &   my_particle,
                                                                             const  double                    my_radius,
                                                                             const  unsigned int              my_number_points_theta,
                                                                             const  unsigned int              my_number_points_phi,
                                                                             const  unsigned int              my_init_size)
: ParticleObserver<FieldModel>(my_particle,my_init_size)
, radius(my_radius)
, number_points_theta(my_number_points_theta)
, number_points_phi(my_number_points_phi)
{
  // Set the sizes of the containers that will store the Liénard-Wiechert fields.
}

template <class FieldModel>
inline
void
ParticleObserverLienardWiechert::operator(const arma::colvec::fixed<8> & x,
                                                double                   t)
{
  // We compute the usual stuff.
  ParticleObserver<FieldModel>::operator(x,t);

  // We compute the Liénard-Wichert fields.
}

template <class FieldModel>
inline
void
ParticleObserverLienardWiechert::WriteAllData(hid_t group_id)
{
  // We write the usual stuff.
  ParticleObserver<FieldModel>::WriteAllData(group_id);

  // We write the Liénard-Wiechert fields.
  WriteLienardWiechertFields(group_id);
}

template <class FieldModel>
inline
void
ParticleObserverLienardWiechert::WriteLienardWiechertFields(hid_t group_id)
{
  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Properties of the datasets.
  // hsize_t size_dataspace[2] = ???;
  // hsize_t chunk_size[2]     = {50,3};
  // SetHDF5Properties(dataspace_id, plist_id, 2, size_dataspace, chunk_size, 5);
}

} // namespace mellotron

#endif // PARTICLE_OBSERVER_MEAT_HPP

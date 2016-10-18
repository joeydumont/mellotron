template <class FieldModel>
inline
ParticleObserver<FieldModel>::ParticleObserver(Particle<FieldModel> & my_particle)
: particle(my_particle)
{}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::operator() (const arma::colvec::fixed<8> &x,
                                                double                  t)
{
  // Resize the matrices.
  const int n_cols = position.n_cols;
  position.resize(3,n_cols+1);
  momentum.resize(3,n_cols+1);
  electric_field.resize(3,n_cols+1);
  magnetic_field.resize(3,n_cols+1);

  // Push the data to the appropriate containers.
  position.col(n_cols) = x.subvec(1,3);
  momentum.col(n_cols) = x.subvec(5,7);
  gamma.push_back(x(4));
  times.push_back(x[0]);

  // Compute the electromagnetic field and store it.
  particle.ComputeFieldTensor(x(0),x(1),x(2),x(3));
  electric_field.col(n_cols) = particle.GetElectricField();
  magnetic_field.col(n_cols) = particle.GetMagneticField();

  // Compute chi.
  double pdotE         = arma::dot(momentum.col(n_cols),electric_field.col(n_cols));
  arma::colvec lorentz = x[4]*electric_field.col(n_cols) + arma::cross(momentum.col(n_cols),magnetic_field.col(n_cols));
  double chi_sq        = std::pow(arma::norm(lorentz,2),2)-std::pow(pdotE,2);
  chi.push_back(std::sqrt(chi_sq));
  std::cout << "Lorentz norm: " << std::pow(arma::norm(lorentz,2),2) << std:: endl;
  std::cout << "p dot E: " << pdotE << std::endl;
  std::cout << "chi_sq: " << chi_sq << std::endl;
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::OutputData()
{
  GenerateXDMF();
  GenerateHDF5();
}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::GenerateXDMF()
{

}

template <class FieldModel>
inline
void
ParticleObserver<FieldModel>::GenerateHDF5()
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

  // IDs used for the datasets.
  hid_t dataspace_id, plist_id, dataset_id;
  herr_t status;

  // Creation of the temporal dataset.
  hsize_t size_dataspace[1];size_dataspace[0] = times.size();
  dataspace_id = H5Screate_simple(1, size_dataspace, NULL);
  plist_id     = H5Pcreate(H5P_DATASET_CREATE);
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

  H5Dclose(dataset_id);
  H5Sclose(dataspace_id);

  // Creation of the position dataset.
  hsize_t size_dataspace_position[2];
  size_dataspace_position[0] = position.n_cols;
  size_dataspace_position[1] = position.n_rows;
  dataspace_id               = H5Screate_simple(2, size_dataspace_position, NULL);
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

  // We close the HDF5 objects.
  H5Pclose(plist_id);
  H5Sclose(dataspace_id);
  H5Gclose(group_id);
  H5Fclose(file_id);
}